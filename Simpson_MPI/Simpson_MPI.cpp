#include "SimpsonRule.h"
#include <chrono>
#include <iomanip>

double f(double x) {
	return 2 * x + 3 * sqrt(x);
}

int main() {
	int rank, comm_size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	SimpsonRule simpson_helper;
	simpson_helper.init(0, 5000, 3e9, rank, comm_size);

	std::chrono::steady_clock::time_point start, end;
	
	start = std::chrono::high_resolution_clock::now();
	double res = simpson_helper.parallel_simpson(f);
	end = std::chrono::high_resolution_clock::now();
	if (rank == 0) {
		double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
		time_taken *= 1e-9;
		std::cout << "result: " << res << " time : " << std::fixed << time_taken << std::setprecision(9) << " sec" << std::endl;
	}

	MPI_Finalize();
	return 0;

}
