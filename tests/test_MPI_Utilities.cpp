#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "MPI_Utilities.hpp"

using namespace DaMaSCUS_SUN;

namespace
{
void Check(
	bool condition,
	const std::string& message,
	int rank,
	int& failures)
{
	if(condition)
		return;
	std::cerr << "[mpi-utilities rank " << rank << "] "
	          << message << std::endl;
	failures++;
}

void Compute_With_Progress(
	MPIWorkQueue& queue,
	std::chrono::milliseconds duration)
{
	const auto deadline = std::chrono::steady_clock::now() + duration;
	while(std::chrono::steady_clock::now() < deadline)
	{
		queue.Progress();
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}
}
}

int main(int argc, char* argv[])
{
	if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
		return 1;

	int rank = 0;
	int processes = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &processes);
	int failures = 0;
	Check(
	    processes == 4,
	    "test requires exactly four MPI ranks",
	    rank,
	    failures);

	int preflight_failure = failures == 0 ? 0 : 1;
	MPI_Allreduce(
	    MPI_IN_PLACE,
	    &preflight_failure,
	    1,
	    MPI_INT,
	    MPI_MAX,
	    MPI_COMM_WORLD);
	if(preflight_failure != 0)
	{
		MPI_Finalize();
		return 1;
	}

	const std::string all_empty =
	    Gather_MPI_Text_To_Root(std::string());
	Check(
	    all_empty.empty(),
	    "all-empty gather did not return an empty string",
	    rank,
	    failures);
	// This barrier catches an implementation that returns from the zero-byte
	// gather on root while leaving non-root ranks inside the collective.
	MPI_Barrier(MPI_COMM_WORLD);

	std::string mixed_local;
	if(rank == 0)
		mixed_local = "rank-0\n";
	else if(rank == 2)
		mixed_local = "rank-2\n";
	const std::string mixed =
	    Gather_MPI_Text_To_Root(mixed_local);
	if(rank == 0)
	{
		Check(
		    mixed == "rank-0\nrank-2\n",
		    "mixed empty/non-empty gather lost rank order or content",
		    rank,
		    failures);
	}
	else
	{
		Check(
		    mixed.empty(),
		    "non-root rank received gathered text",
		    rank,
		    failures);
	}

	const std::string every_rank_local =
	    std::to_string(rank) + ",";
	const std::string nonzero_root =
	    Gather_MPI_Text_To_Root(every_rank_local, 2);
	if(rank == 2)
	{
		Check(
		    nonzero_root == "0,1,2,3,",
		    "nonzero-root gather lost rank order or content",
		    rank,
		    failures);
	}
	else
	{
		Check(
		    nonzero_root.empty(),
		    "rank other than the selected root received gathered text",
		    rank,
		    failures);
	}

	// A numerical failure is recorded and replaced; it must not stop the queue.
	MPI_Barrier(MPI_COMM_WORLD);
	MPIWorkQueue replacement_queue(
	    4,
	    8,
	    1.0,
	    MPI_COMM_WORLD);
	// Reserve one task on every rank before completing any of them, so rank 0
	// always supplies the intended failure regardless of process scheduling.
	const MPIWorkClaimResult initial_replacement_claim =
	    replacement_queue.TryClaim();
	Check(
	    initial_replacement_claim == MPIWorkClaimResult::Claimed,
	    "rank did not receive its initial replacement-test task",
	    rank,
	    failures);
	MPI_Barrier(MPI_COMM_WORLD);
	if(initial_replacement_claim == MPIWorkClaimResult::Claimed)
	{
		MPIWorkOutcome outcome;
		outcome.numerical_failure = rank == 0;
		outcome.accepted_sample = rank != 0;
		replacement_queue.Complete(outcome);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	while(true)
	{
		const MPIWorkClaimResult claim =
		    replacement_queue.TryClaim();
		if(claim == MPIWorkClaimResult::Stop)
			break;
		if(claim == MPIWorkClaimResult::Wait)
		{
			std::this_thread::sleep_for(
			    std::chrono::milliseconds(1));
			continue;
		}

		MPIWorkOutcome outcome;
		outcome.accepted_sample = true;
		replacement_queue.Complete(outcome);
	}
	const MPIWorkQueueState replacement_state =
	    replacement_queue.Finalize();
	Check(
	    replacement_state.work_claims == 5
	        && replacement_state.completed_trajectories == 5
	        && replacement_state.accepted_samples == 4
	        && replacement_state.numerical_failures == 1
	        && replacement_state.stop_reason
	           == MPIWorkStopReason::None,
	    "numerical failure was not recorded and replaced",
	    rank,
	    failures);

	// Synthetic imbalance regression: the slow rank claims one long task while
	// faster ranks immediately recycle their completed slots. Making a non-root
	// rank slow also verifies that its final RMA completion can reach rank 0
	// after faster ranks have entered collective queue finalization.
	MPI_Barrier(MPI_COMM_WORLD);
	MPIWorkQueue work_queue(
	    8,
	    64,
	    1.0,
	    MPI_COMM_WORLD);
	uint64_t local_completed = 0;
	auto complete_imbalanced_task = [&]()
	{
		Compute_With_Progress(
		    work_queue,
		    std::chrono::milliseconds(rank == 3 ? 80 : 5));
		MPIWorkOutcome outcome;
		outcome.accepted_sample = true;
		work_queue.Complete(outcome);
		local_completed++;
	};
	const MPIWorkClaimResult initial_imbalanced_claim =
	    work_queue.TryClaim();
	Check(
	    initial_imbalanced_claim == MPIWorkClaimResult::Claimed,
	    "rank did not receive its initial imbalance-test task",
	    rank,
	    failures);
	MPI_Barrier(MPI_COMM_WORLD);
	if(initial_imbalanced_claim == MPIWorkClaimResult::Claimed)
		complete_imbalanced_task();
	while(true)
	{
		const MPIWorkClaimResult claim =
		    work_queue.TryClaim();
		if(claim == MPIWorkClaimResult::Stop)
			break;
		if(claim == MPIWorkClaimResult::Wait)
		{
			std::this_thread::sleep_for(
			    std::chrono::milliseconds(1));
			continue;
		}

		complete_imbalanced_task();
	}
	const MPIWorkQueueState queue_state =
	    work_queue.Finalize();
	Check(
	    queue_state.work_claims == 8
	        && queue_state.completed_trajectories == 8
	        && queue_state.accepted_samples == 8
	        && queue_state.in_flight == 0,
	    "dynamic queue did not close exactly at the requested target",
	    rank,
	    failures);
	Check(
	    queue_state.peak_in_flight == 4,
	    "dynamic queue did not keep all four ranks busy initially",
	    rank,
	    failures);
	Check(
	    queue_state.stop_reason == MPIWorkStopReason::None,
	    "successful dynamic queue reported an early-stop reason",
	    rank,
	    failures);

	std::vector<uint64_t> completed_by_rank(4, 0);
	MPI_Allgather(
	    &local_completed,
	    1,
	    MPI_UINT64_T,
	    completed_by_rank.data(),
	    1,
	    MPI_UINT64_T,
	    MPI_COMM_WORLD);
	Check(
	    completed_by_rank[3] == 1,
	    "slow rank received more work while faster ranks were available",
	    rank,
	    failures);
	Check(
	    *std::max_element(
	        completed_by_rank.begin(),
	        completed_by_rank.begin() + 3) > 1,
	    "fast ranks did not recycle work dynamically",
	    rank,
	    failures);

	// Rank 0 retains a long-running task while every other rank starts and
	// completes two tasks. During this interval root calls only Progress():
	// polling a notification with MPI_Test/Iprobe here would itself advance
	// RMA and conceal a broken Progress() implementation. Workers record their
	// elapsed time before the later gather, which also lets a failing test drain
	// safely once root's bounded computation interval has ended.
	MPIWorkQueue slow_root_queue(7, 7, 1.0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		Check(
		    slow_root_queue.TryClaim() == MPIWorkClaimResult::Claimed,
		    "root could not reserve its long-running task",
		    rank,
		    failures);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	const auto slow_root_start = std::chrono::steady_clock::now();
	double worker_completed_after_sec = 0.0;
	if(rank == 0)
	{
		Compute_With_Progress(
		    slow_root_queue,
		    std::chrono::milliseconds(600));
	}
	else
	{
		for(int task = 0; task < 2; task++)
		{
			const MPIWorkClaimResult claim = slow_root_queue.TryClaim();
			Check(
			    claim == MPIWorkClaimResult::Claimed,
			    "worker could not recycle work while root was computing",
			    rank,
			    failures);
			if(claim != MPIWorkClaimResult::Claimed)
				break;
			MPIWorkOutcome outcome;
			outcome.accepted_sample = true;
			slow_root_queue.Complete(outcome);
		}
		worker_completed_after_sec = std::chrono::duration<double>(
		    std::chrono::steady_clock::now() - slow_root_start).count();
	}
	std::vector<double> worker_completion_times(4, 0.0);
	MPI_Gather(
	    &worker_completed_after_sec,
	    1,
	    MPI_DOUBLE,
	    worker_completion_times.data(),
	    1,
	    MPI_DOUBLE,
	    0,
	    MPI_COMM_WORLD);
	if(rank == 0)
	{
		for(int worker = 1; worker < processes; worker++)
		{
			Check(
			    worker_completion_times[static_cast<size_t>(worker)] < 0.4,
			    "worker " + std::to_string(worker)
			        + " waited for root's long computation before completing work",
			    rank,
			    failures);
		}
		const MPIWorkQueueState before_root_completion =
		    slow_root_queue.ReadState();
		Check(
		    before_root_completion.work_claims == 7
		        && before_root_completion.completed_trajectories == 6
		        && before_root_completion.accepted_samples == 6
		        && before_root_completion.in_flight == 1,
		    "workers did not complete multiple tasks while root retained its claim",
		    rank,
		    failures);
		MPIWorkOutcome outcome;
		outcome.accepted_sample = true;
		slow_root_queue.Complete(outcome);
	}
	const MPIWorkQueueState slow_root_final = slow_root_queue.Finalize();
	Check(
	    slow_root_final.accepted_samples == 7
	        && slow_root_final.completed_trajectories == 7
	        && slow_root_final.in_flight == 0,
	    "slow-root queue did not drain after root completed",
	    rank,
	    failures);

	// Exact target semantics intentionally limit parallelism when fewer target
	// slots remain than MPI ranks. A target of one must never issue a second
	// claim while the first trajectory is still running.
	MPIWorkQueue one_sample_queue(1, 8, 1.0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		Check(
		    one_sample_queue.TryClaim() == MPIWorkClaimResult::Claimed,
		    "root could not reserve the only target slot",
		    rank,
		    failures);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank != 0)
	{
		MPIWorkQueueState observed;
		Check(
		    one_sample_queue.TryClaim(&observed) == MPIWorkClaimResult::Wait
		        && observed.work_claims == 1
		        && observed.in_flight == 1
		        && observed.accepted_samples == 0,
		    "target-one queue issued excess work instead of waiting",
		    rank,
		    failures);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
	{
		MPIWorkOutcome outcome;
		outcome.accepted_sample = true;
		one_sample_queue.Complete(outcome);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	Check(
	    one_sample_queue.TryClaim() == MPIWorkClaimResult::Stop,
	    "target-one queue did not stop after its accepted sample",
	    rank,
	    failures);
	const MPIWorkQueueState one_sample_final = one_sample_queue.Finalize();
	Check(
	    one_sample_final.work_claims == 1
	        && one_sample_final.completed_trajectories == 1
	        && one_sample_final.accepted_samples == 1
	        && one_sample_final.peak_in_flight == 1
	        && one_sample_final.in_flight == 0,
	    "target-one queue failed exact-target accounting",
	    rank,
	    failures);

	int global_failures = 0;
	MPI_Allreduce(
	    &failures,
	    &global_failures,
	    1,
	    MPI_INT,
	    MPI_SUM,
	    MPI_COMM_WORLD);
	MPI_Finalize();
	return global_failures == 0 ? 0 : 1;
}
