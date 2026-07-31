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

	// Synthetic imbalance regression: the slow rank claims one long task while
	// faster ranks immediately recycle their completed slots. Making a non-root
	// rank slow also verifies that its final RMA completion can reach rank 0
	// after faster ranks have entered collective queue finalization.
	MPI_Barrier(MPI_COMM_WORLD);
	MPIWorkQueue work_queue(
	    8,
	    64,
	    false,
	    1.0,
	    1.0,
	    MPI_COMM_WORLD);
	uint64_t local_completed = 0;
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

		std::this_thread::sleep_for(
		    std::chrono::milliseconds(rank == 3 ? 80 : 5));
		MPIWorkOutcome outcome;
		outcome.accepted_sample = true;
		work_queue.Complete(outcome);
		local_completed++;
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
