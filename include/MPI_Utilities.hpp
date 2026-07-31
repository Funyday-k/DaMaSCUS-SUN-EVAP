#ifndef DAMASCUS_SUN_MPI_UTILITIES_HPP
#define DAMASCUS_SUN_MPI_UTILITIES_HPP

#include <mpi.h>

#include <string>

namespace DaMaSCUS_SUN
{

// Collectively concatenate rank-local text in rank order. The concatenated
// text is returned only on root; every other rank returns an empty string.
//
// The all-empty case deliberately skips MPI_Gatherv. Some MPI implementations
// do not reliably return on every rank when a zero-count Gatherv is called with
// null buffers. For mixed empty/non-empty inputs, valid dummy buffers are used
// on zero-length ranks.
std::string Gather_MPI_Text_To_Root(
	const std::string& local_text,
	int root = 0,
	MPI_Comm communicator = MPI_COMM_WORLD);

}	// namespace DaMaSCUS_SUN

#endif
