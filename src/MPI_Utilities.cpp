#include "MPI_Utilities.hpp"

#include <algorithm>
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

}	// namespace DaMaSCUS_SUN
