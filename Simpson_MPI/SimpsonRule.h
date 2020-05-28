#pragma once
#include<iostream>
#include <omp.h>
#include<mpi.h>
class SimpsonRule {
	long long n;
	double xmin, xmax, step;
	int _rank, _comm_size;
public:
	double parallel_simpson(double (*func)(double x));
	void init(double xmin, double xmax, long long n, int rank, int comm_size);
};
