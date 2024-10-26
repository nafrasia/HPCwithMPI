#ifndef MPISTUFF_H_HEADER
#define MPISTUFF_H_HEADER

#include "mpi.h"

class MPI_Stuff
{
	public:
		int nID, nprocs;

		MPI_Stuff(int & argc, char** &argv)
		{
			MPI_Init(&argc, &argv);
			MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
			MPI_Comm_rank(MPI_COMM_WORLD, &nID);
		};
		~MPI_Stuff()
		{
			MPI_Finalize();
		};
};

#endif
