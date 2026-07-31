#ifndef DAMASCUS_SUN_MPI_UTILITIES_HPP
#define DAMASCUS_SUN_MPI_UTILITIES_HPP

#include <mpi.h>

#include <string>

namespace DaMaSCUS_SUN
{

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

}	// namespace DaMaSCUS_SUN

#endif
