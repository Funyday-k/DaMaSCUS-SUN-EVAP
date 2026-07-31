#include <mpi.h>

#include <iostream>
#include <string>

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
