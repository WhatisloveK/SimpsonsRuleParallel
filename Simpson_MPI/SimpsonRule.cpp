#include "SimpsonRule.h"


void SimpsonRule::init(double xmin, double xmax, long long n, int rank, int comm_size) {
	this->xmax = xmax;
	this->xmin = xmin;
	this->n = n;
	this->step = (xmax - xmin) / n;
	this->_rank = rank;
	this->_comm_size = comm_size;
}

double SimpsonRule::parallel_simpson(double (*func)(double x)) {
	double simpson_integral = 0;
	double sum2 = 0;
	double sum4 = 0;
	double sum = 0;

	long long length = n / _comm_size;

	long long left = 1 + _rank * length;
	long long right = __min(n, left + length - 1);


	for (long long i = left; i <= right; i += 2)
	{
		sum4 += func(xmin + step * i);
		sum2 += func(xmin + step * (i + 1));
	}
	

	double SUM4 = 0;
	double SUM2 = 0;

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(&sum4, &SUM4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sum2, &SUM2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (_rank == 0) {
		sum = func(xmin) + 4 * SUM4 + 2 * SUM2 - func(xmax);
		simpson_integral = (step / 3) * sum;
		return simpson_integral;
	}

	return 0;
}