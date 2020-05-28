#include<iostream>
#include<mpi.h>

int main() {
	int rank, comm_size, t;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	if (rank == 0) {
		std::cin >> t;
	}


	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;

}
