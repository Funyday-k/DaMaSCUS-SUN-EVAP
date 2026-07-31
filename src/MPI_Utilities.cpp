#include "MPI_Utilities.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace DaMaSCUS_SUN
{

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
	uint64_t total_bytes_u64 = 0;
	if(MPI_Allreduce(
	       &local_bytes_u64,
	       &total_bytes_u64,
	       1,
	       MPI_UINT64_T,
	       MPI_SUM,
	       communicator)
	   != MPI_SUCCESS)
	{
		throw std::runtime_error(
		    "Gather_MPI_Text_To_Root(): failed to total rank-local text sizes.");
	}
	if(total_bytes_u64
	   > static_cast<uint64_t>(std::numeric_limits<int>::max()))
	{
		throw std::overflow_error(
		    "Gather_MPI_Text_To_Root(): gathered text exceeds MPI's int count range.");
	}

	// Avoid the zero-count/null-buffer Gatherv path entirely. Intel MPI 2021.6
	// can return on root while leaving non-root ranks spinning in this case.
	if(total_bytes_u64 == 0)
		return std::string();

	const int local_bytes = static_cast<int>(local_bytes_u64);
	std::vector<int> byte_counts(
	    mpi_rank == root ? static_cast<size_t>(mpi_processes) : 0,
	    0);
	int ignored_count = 0;
	if(MPI_Gather(
	       &local_bytes,
	       1,
	       MPI_INT,
	       mpi_rank == root ? byte_counts.data() : &ignored_count,
	       1,
	       MPI_INT,
	       root,
	       communicator)
	   != MPI_SUCCESS)
	{
		throw std::runtime_error(
		    "Gather_MPI_Text_To_Root(): failed to gather rank-local text sizes.");
	}

	std::vector<int> displacements(
	    mpi_rank == root ? static_cast<size_t>(mpi_processes) : 0,
	    0);
	if(mpi_rank == root)
	{
		for(int rank = 1; rank < mpi_processes; rank++)
			displacements[static_cast<size_t>(rank)] =
			    displacements[static_cast<size_t>(rank - 1)]
			    + byte_counts[static_cast<size_t>(rank - 1)];
	}

	std::vector<char> gathered_text(
	    mpi_rank == root ? static_cast<size_t>(total_bytes_u64) : 0);
	char ignored_byte = '\0';
	const char* send_buffer =
	    local_bytes > 0 ? local_text.data() : &ignored_byte;
	char* receive_buffer =
	    mpi_rank == root ? gathered_text.data() : &ignored_byte;
	int* receive_counts =
	    mpi_rank == root ? byte_counts.data() : &ignored_count;
	int* receive_displacements =
	    mpi_rank == root ? displacements.data() : &ignored_count;

	if(MPI_Gatherv(
	       send_buffer,
	       local_bytes,
	       MPI_CHAR,
	       receive_buffer,
	       receive_counts,
	       receive_displacements,
	       MPI_CHAR,
	       root,
	       communicator)
	   != MPI_SUCCESS)
	{
		throw std::runtime_error(
		    "Gather_MPI_Text_To_Root(): failed to gather rank-local text.");
	}

	if(mpi_rank != root)
		return std::string();
	return std::string(gathered_text.begin(), gathered_text.end());
}

}	// namespace DaMaSCUS_SUN
