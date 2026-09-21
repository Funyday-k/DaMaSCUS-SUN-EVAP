#ifndef DAMASCUS_SUN_MPI_UTILITIES_HPP
#define DAMASCUS_SUN_MPI_UTILITIES_HPP

#include <mpi.h>

#include <cstdint>
#include <chrono>
#include <string>
#include <vector>

namespace DaMaSCUS_SUN
{

// Sum variable-length radial histograms on every rank. Missing outer bins
// contribute zero; the result has the largest rank-local extent.
void Allreduce_MPI_Histogram(
	std::vector<double>& histogram,
	MPI_Comm communicator = MPI_COMM_WORLD);

// Collectively concatenate rank-local byte buffers in rank order. The
// concatenated bytes are returned only on root; every other rank returns an
// empty vector. Empty buffers never enter a zero-count data operation.
std::vector<char> Gather_MPI_Bytes_To_Root(
	const void* local_data,
	uint64_t local_bytes,
	int root = 0,
	MPI_Comm communicator = MPI_COMM_WORLD);

// Collectively concatenate rank-local text in rank order. The concatenated
// text is returned only on root; every other rank returns an empty string.
//
// This deliberately avoids MPI_Gatherv. Intel MPI 2021.6 can leave a rank
// spinning in variable-count gathers when one or more inputs are empty, even
// when valid dummy buffers are supplied. Sizes are exchanged with a fixed-count
// collective, and root receives data only from ranks with non-empty text.
std::string Gather_MPI_Text_To_Root(
	const std::string& local_text,
	int root = 0,
	MPI_Comm communicator = MPI_COMM_WORLD);

enum class MPIWorkStopReason : uint64_t
{
	None = 0,
	MaxTrajectoriesReached = 1,
	InitialShiftFailureFractionExceeded = 2,
	RejectedTrajectory = 3
};

struct MPIWorkOutcome
{
	bool accepted_sample = false;
	bool initial_shift_failure = false;
	bool numerical_failure = false;
	bool computational_truncation = false;
	// Production stops issuing work after a failed history; in-flight work drains.
	bool reject_run = false;
};

// Public diagnostic snapshot of the dynamic trajectory queue. Claims are
// globally bounded so accepted_samples + in_flight never exceeds the requested
// target; this preserves the exact target without static rank-sized batches.
struct MPIWorkQueueState
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
	MPIWorkStopReason stop_reason = MPIWorkStopReason::None;
};

enum class MPIWorkClaimResult
{
	Claimed,
	Wait,
	Stop
};

// MPI RMA trajectory queue used by every rank. Each completion immediately
// releases capacity for another claim, so a fast rank can keep working without
// waiting for the slowest trajectory in a collective batch.
class MPIWorkQueue
{
  public:
	MPIWorkQueue(
		uint64_t target_samples,
		uint64_t maximum_trajectories,
		double initial_shift_abort_fraction,
		MPI_Comm communicator = MPI_COMM_WORLD);
	~MPIWorkQueue();

	MPIWorkQueue(const MPIWorkQueue&) = delete;
	MPIWorkQueue& operator=(const MPIWorkQueue&) = delete;

	MPIWorkClaimResult TryClaim(MPIWorkQueueState* observed_state = nullptr);
	MPIWorkQueueState Complete(const MPIWorkOutcome& outcome);
	MPIWorkQueueState ReadState() const;

	// Call periodically from the main thread while computing a trajectory.
	// Passive-target RMA may need rank 0 to enter MPI to serve other ranks.
	// This never takes the queue lock, waits for a peer, or changes work state.
	void Progress();

	// Collective. Every rank must stop claiming/completing work before calling.
	MPIWorkQueueState Finalize();

  private:
	struct WireState;
	static MPIWorkQueueState PublicState(const WireState& state);
	WireState ReadWireStateLocked() const;
	void WriteWireStateLocked(const WireState& state) const;
	void LockRoot() const;
	void UnlockRoot() const;

	MPI_Comm communicator_;
	int mpi_rank_ = 0;
	int mpi_processes_ = 0;
	double initial_shift_abort_fraction_ = 0.0;
	mutable MPI_Win window_ = MPI_WIN_NULL;
	void* root_window_memory_ = nullptr;
	bool finalized_ = false;
	std::chrono::steady_clock::time_point next_progress_poll_{};
};

}	// namespace DaMaSCUS_SUN

#endif
