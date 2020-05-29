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
	double sub_sum2 = 0;
	double sub_sum4 = 0;
	double sum = 0;

	long long left; 
	long long right;

	if (_rank == 0) {
		long long length;
		length = n / _comm_size;
		
		for (int i = 1; i < _comm_size; i++) {
			left = i * length;
			right = __min(n, left + length);
			MPI_Send(&left,
				1, MPI_DOUBLE,
				i, i,
				MPI_COMM_WORLD);
			MPI_Send(&right,
				1, MPI_DOUBLE,
				i, i,
				MPI_COMM_WORLD);
		}

		left = 0;
		right = length;
		for (long long i = left; i <= right; i += 2)
		{
			sub_sum4 += func(xmin + step * i);
			sub_sum2 += func(xmin + step * (i + 1));
		}



		double SUB_SUM4, SUB_SUM2, SUM4, SUM2;

		SUM4 = sub_sum4;
		SUM2 = sub_sum2;
		for (int i = 1; i < _comm_size; i++) {
			MPI_Recv(&SUB_SUM2, 1, MPI_DOUBLE,
				i, i,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);
			MPI_Recv(&SUB_SUM4, 1, MPI_DOUBLE,
				i, i,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);
			SUM4 += SUB_SUM4;
			SUM2 += SUB_SUM2;
		}

		sum = func(xmin) + 4 * SUM4 + 2 * SUM2 - func(xmax);
		simpson_integral = (step / 3) * sum;
		return simpson_integral;
	}
	else {

		MPI_Recv(&left, 1, MPI_DOUBLE,
			0, _rank,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		MPI_Recv(&right, 1, MPI_DOUBLE,
			0, _rank,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		for (long long i = left; i <= right; i += 2)
		{
			sub_sum4 += func(xmin + step * i);
			sub_sum2 += func(xmin + step * (i + 1));
		}
		MPI_Send(&sub_sum2,
			1, MPI_DOUBLE,
			0, _rank,
			MPI_COMM_WORLD);
		MPI_Send(&sub_sum4,
			1, MPI_DOUBLE,
			0, _rank,
			MPI_COMM_WORLD);
	}

	return 0;
}