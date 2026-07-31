#include "MPI_Utilities.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace DaMaSCUS_SUN
{
namespace
{
void Check_MPI_Result(int result, const char* operation)
{
	if(result != MPI_SUCCESS)
		throw std::runtime_error(
		    std::string("MPI work queue failed during ") + operation + ".");
}
}

struct MPIWorkQueue::WireState
{
	uint64_t target_samples = 0;
	uint64_t maximum_trajectories = 0;
	uint64_t work_claims = 0;
	uint64_t completed_trajectories = 0;
	uint64_t accepted_samples = 0;
	uint64_t in_flight = 0;
	uint64_t peak_in_flight = 0;
	uint64_t initial_shift_failures = 0;
	uint64_t numerical_failures = 0;
	uint64_t computational_truncations = 0;
	uint64_t stop_reason = 0;
};

std::string Gather_MPI_Text_To_Root(
	const std::string& local_text,
	int root,
	MPI_Comm communicator)
{
	int mpi_rank = 0;
	int mpi_processes = 0;
	if(MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS
	   || MPI_Comm_size(communicator, &mpi_processes) != MPI_SUCCESS)
	{
		throw std::runtime_error(
		    "Gather_MPI_Text_To_Root(): failed to inspect the MPI communicator.");
	}
	if(root < 0 || root >= mpi_processes)
		throw std::invalid_argument(
		    "Gather_MPI_Text_To_Root(): root rank is outside the MPI communicator.");

	const uint64_t local_bytes_u64 =
	    static_cast<uint64_t>(local_text.size());
	std::vector<uint64_t> byte_counts_u64(
	    static_cast<size_t>(mpi_processes),
	    0);
	if(MPI_Allgather(
	       &local_bytes_u64,
	       1,
	       MPI_UINT64_T,
	       byte_counts_u64.data(),
	       1,
	       MPI_UINT64_T,
	       communicator)
	   != MPI_SUCCESS)
	{
		throw std::runtime_error(
		    "Gather_MPI_Text_To_Root(): failed to exchange rank-local text sizes.");
	}

	uint64_t total_bytes_u64 = 0;
	for(uint64_t rank_bytes : byte_counts_u64)
	{
		if(rank_bytes
		   > static_cast<uint64_t>(std::numeric_limits<int>::max())
		         - total_bytes_u64)
		{
			throw std::overflow_error(
			    "Gather_MPI_Text_To_Root(): gathered text exceeds MPI's int count range.");
		}
		total_bytes_u64 += rank_bytes;
	}

	// Avoid any zero-count data operation in the all-empty case.
	if(total_bytes_u64 == 0)
		return std::string();

	const int local_bytes = static_cast<int>(local_bytes_u64);
	std::vector<int> byte_counts(
	    static_cast<size_t>(mpi_processes),
	    0);
	std::vector<int> displacements(
	    static_cast<size_t>(mpi_processes),
	    0);
	for(int rank = 0; rank < mpi_processes; rank++)
	{
		byte_counts[static_cast<size_t>(rank)] =
		    static_cast<int>(byte_counts_u64[static_cast<size_t>(rank)]);
		if(rank > 0)
			displacements[static_cast<size_t>(rank)] =
			    displacements[static_cast<size_t>(rank - 1)]
			    + byte_counts[static_cast<size_t>(rank - 1)];
	}

	std::vector<char> gathered_text(
	    mpi_rank == root ? static_cast<size_t>(total_bytes_u64) : 0);
	constexpr int text_gather_tag = 7319;
	if(mpi_rank == root)
	{
		if(local_bytes > 0)
		{
			std::copy(
			    local_text.begin(),
			    local_text.end(),
			    gathered_text.begin()
			        + displacements[static_cast<size_t>(root)]);
		}
		for(int rank = 0; rank < mpi_processes; rank++)
		{
			const int rank_bytes =
			    byte_counts[static_cast<size_t>(rank)];
			if(rank == root || rank_bytes == 0)
				continue;
			if(MPI_Recv(
			       gathered_text.data()
			           + displacements[static_cast<size_t>(rank)],
			       rank_bytes,
			       MPI_CHAR,
			       rank,
			       text_gather_tag,
			       communicator,
			       MPI_STATUS_IGNORE)
			   != MPI_SUCCESS)
			{
				throw std::runtime_error(
				    "Gather_MPI_Text_To_Root(): failed to receive rank-local text.");
			}
		}
	}
	else if(local_bytes > 0
	        && MPI_Send(
	               local_text.data(),
	               local_bytes,
	               MPI_CHAR,
	               root,
	               text_gather_tag,
	               communicator)
	               != MPI_SUCCESS)
	{
		throw std::runtime_error(
		    "Gather_MPI_Text_To_Root(): failed to send rank-local text.");
	}

	// Keep later collectives from overtaking ranks that were still blocked in
	// a send while root received another rank's text.
	if(MPI_Barrier(communicator) != MPI_SUCCESS)
		throw std::runtime_error(
		    "Gather_MPI_Text_To_Root(): failed to synchronize ranks after gathering text.");

	if(mpi_rank != root)
		return std::string();
	return std::string(gathered_text.begin(), gathered_text.end());
}

MPIWorkQueue::MPIWorkQueue(
	uint64_t target_samples,
	uint64_t maximum_trajectories,
	bool abort_on_invalid_fraction,
	double initial_shift_abort_fraction,
	double invalid_abort_fraction,
	MPI_Comm communicator)
: communicator_(communicator),
  abort_on_invalid_fraction_(abort_on_invalid_fraction),
  initial_shift_abort_fraction_(initial_shift_abort_fraction),
  invalid_abort_fraction_(invalid_abort_fraction)
{
	if(target_samples == 0)
		throw std::invalid_argument(
		    "MPIWorkQueue(): target_samples must be positive.");
	if(maximum_trajectories == 0)
		throw std::invalid_argument(
		    "MPIWorkQueue(): maximum_trajectories must be positive.");
	if(!std::isfinite(initial_shift_abort_fraction_)
	   || initial_shift_abort_fraction_ < 0.0
	   || !std::isfinite(invalid_abort_fraction_)
	   || invalid_abort_fraction_ < 0.0)
	{
		throw std::invalid_argument(
		    "MPIWorkQueue(): abort fractions must be finite and non-negative.");
	}

	Check_MPI_Result(
	    MPI_Comm_rank(communicator_, &mpi_rank_),
	    "communicator-rank inspection");
	Check_MPI_Result(
	    MPI_Comm_size(communicator_, &mpi_processes_),
	    "communicator-size inspection");
	if(mpi_processes_ <= 0)
		throw std::runtime_error(
		    "MPIWorkQueue(): communicator contains no ranks.");

	const MPI_Aint local_window_bytes =
	    mpi_rank_ == 0
	    ? static_cast<MPI_Aint>(sizeof(WireState))
	    : static_cast<MPI_Aint>(0);
	Check_MPI_Result(
	    MPI_Win_allocate(
	        local_window_bytes,
	        1,
	        MPI_INFO_NULL,
	        communicator_,
	        &root_window_memory_,
	        &window_),
	    "window allocation");

	if(mpi_rank_ == 0)
	{
		WireState initial;
		initial.target_samples = target_samples;
		initial.maximum_trajectories = maximum_trajectories;
		LockRoot();
		try
		{
			WriteWireStateLocked(initial);
			UnlockRoot();
		}
		catch(...)
		{
			MPI_Win_unlock(0, window_);
			throw;
		}
	}
	Check_MPI_Result(
	    MPI_Barrier(communicator_),
	    "queue initialization barrier");
}

MPIWorkQueue::~MPIWorkQueue()
{
	// MPI_Win_free is collective and therefore cannot safely be called from a
	// destructor during exception unwinding. Normal callers finish with
	// Finalize(), which releases the window on every rank.
}

MPIWorkQueueState MPIWorkQueue::PublicState(const WireState& state)
{
	MPIWorkQueueState result;
	result.target_samples = state.target_samples;
	result.maximum_trajectories = state.maximum_trajectories;
	result.work_claims = state.work_claims;
	result.completed_trajectories = state.completed_trajectories;
	result.accepted_samples = state.accepted_samples;
	result.in_flight = state.in_flight;
	result.peak_in_flight = state.peak_in_flight;
	result.initial_shift_failures = state.initial_shift_failures;
	result.numerical_failures = state.numerical_failures;
	result.computational_truncations =
	    state.computational_truncations;
	result.stop_reason =
	    static_cast<MPIWorkStopReason>(state.stop_reason);
	return result;
}

void MPIWorkQueue::LockRoot() const
{
	Check_MPI_Result(
	    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, window_),
	    "root-state lock");
}

void MPIWorkQueue::UnlockRoot() const
{
	Check_MPI_Result(
	    MPI_Win_unlock(0, window_),
	    "root-state unlock");
}

MPIWorkQueue::WireState MPIWorkQueue::ReadWireStateLocked() const
{
	WireState state;
	Check_MPI_Result(
	    MPI_Get(
	        &state,
	        static_cast<int>(sizeof(state)),
	        MPI_BYTE,
	        0,
	        0,
	        static_cast<int>(sizeof(state)),
	        MPI_BYTE,
	        window_),
	    "root-state read");
	Check_MPI_Result(
	    MPI_Win_flush(0, window_),
	    "root-state read completion");
	return state;
}

void MPIWorkQueue::WriteWireStateLocked(const WireState& state) const
{
	Check_MPI_Result(
	    MPI_Put(
	        &state,
	        static_cast<int>(sizeof(state)),
	        MPI_BYTE,
	        0,
	        0,
	        static_cast<int>(sizeof(state)),
	        MPI_BYTE,
	        window_),
	    "root-state write");
}

MPIWorkClaimResult MPIWorkQueue::TryClaim(
	MPIWorkQueueState* observed_state)
{
	if(finalized_)
		throw std::logic_error(
		    "MPIWorkQueue::TryClaim(): queue is already finalized.");

	LockRoot();
	WireState state;
	MPIWorkClaimResult result = MPIWorkClaimResult::Wait;
	try
	{
		state = ReadWireStateLocked();
		if(state.stop_reason
		       != static_cast<uint64_t>(MPIWorkStopReason::None)
		   || state.accepted_samples >= state.target_samples)
		{
			result = MPIWorkClaimResult::Stop;
		}
		else if(state.work_claims >= state.maximum_trajectories)
		{
			if(state.in_flight == 0)
			{
				state.stop_reason = static_cast<uint64_t>(
				    MPIWorkStopReason::MaxTrajectoriesReached);
				WriteWireStateLocked(state);
				result = MPIWorkClaimResult::Stop;
			}
		}
		else if(state.accepted_samples + state.in_flight
		        < state.target_samples)
		{
			state.work_claims++;
			state.in_flight++;
			state.peak_in_flight =
			    std::max(state.peak_in_flight, state.in_flight);
			WriteWireStateLocked(state);
			result = MPIWorkClaimResult::Claimed;
		}
		UnlockRoot();
	}
	catch(...)
	{
		MPI_Win_unlock(0, window_);
		throw;
	}

	if(observed_state != nullptr)
		*observed_state = PublicState(state);
	return result;
}

MPIWorkQueueState MPIWorkQueue::Complete(
	const MPIWorkOutcome& outcome)
{
	if(finalized_)
		throw std::logic_error(
		    "MPIWorkQueue::Complete(): queue is already finalized.");

	LockRoot();
	WireState state;
	try
	{
		state = ReadWireStateLocked();
		if(state.in_flight == 0
		   || state.completed_trajectories >= state.work_claims)
		{
			throw std::logic_error(
			    "MPIWorkQueue::Complete(): no claimed work is in flight.");
		}

		state.in_flight--;
		state.completed_trajectories++;
		if(outcome.accepted_sample)
			state.accepted_samples++;
		if(outcome.initial_shift_failure)
			state.initial_shift_failures++;
		if(outcome.numerical_failure)
			state.numerical_failures++;
		if(outcome.computational_truncation)
			state.computational_truncations++;
		if(state.accepted_samples > state.target_samples)
			throw std::logic_error(
			    "MPIWorkQueue::Complete(): accepted target overshot.");

		if(state.stop_reason
		   == static_cast<uint64_t>(MPIWorkStopReason::None))
		{
			const double attempted =
			    static_cast<double>(state.completed_trajectories);
			const double initial_shift_failure_fraction =
			    static_cast<double>(state.initial_shift_failures)
			    / attempted;
			const double invalid_fraction =
			    static_cast<double>(
			        state.numerical_failures
			        + state.computational_truncations)
			    / attempted;
			if(initial_shift_failure_fraction
			   > initial_shift_abort_fraction_)
			{
				state.stop_reason = static_cast<uint64_t>(
				    MPIWorkStopReason::
				        InitialShiftFailureFractionExceeded);
			}
			else if(abort_on_invalid_fraction_
			        && invalid_fraction
			           > invalid_abort_fraction_)
			{
				state.stop_reason = static_cast<uint64_t>(
				    MPIWorkStopReason::
				        InvalidTrajectoryFractionExceeded);
			}
		}

		WriteWireStateLocked(state);
		UnlockRoot();
	}
	catch(...)
	{
		MPI_Win_unlock(0, window_);
		throw;
	}
	return PublicState(state);
}

MPIWorkQueueState MPIWorkQueue::ReadState() const
{
	if(finalized_)
		throw std::logic_error(
		    "MPIWorkQueue::ReadState(): queue is already finalized.");
	LockRoot();
	WireState state;
	try
	{
		state = ReadWireStateLocked();
		UnlockRoot();
	}
	catch(...)
	{
		MPI_Win_unlock(0, window_);
		throw;
	}
	return PublicState(state);
}

MPIWorkQueueState MPIWorkQueue::Finalize()
{
	if(finalized_)
		throw std::logic_error(
		    "MPIWorkQueue::Finalize(): queue is already finalized.");

	Check_MPI_Result(
	    MPI_Barrier(communicator_),
	    "queue finalization entry barrier");
	const MPIWorkQueueState final_state = ReadState();
	Check_MPI_Result(
	    MPI_Barrier(communicator_),
	    "queue finalization state barrier");
	Check_MPI_Result(
	    MPI_Win_free(&window_),
	    "window release");
	finalized_ = true;
	root_window_memory_ = nullptr;
	return final_state;
}

}	// namespace DaMaSCUS_SUN
