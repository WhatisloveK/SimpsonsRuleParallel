#include "SimpsonRule.h"
#include <omp.h>

void SimpsonRule::set_step(double xmin, double xmax, long long n) {
	this->xmax = xmax;
	this->xmin = xmin;
	this->n = n;
	this->step = (xmax - xmin) / n;
}

double SimpsonRule::parallel_simpson(double (*func)(double x)) {
	double simpson_integral = 0;

#pragma omp parallel for shared(n, xmin, xmax, step, simpson_integral)
	for (long long i = 0; i < n; i++) {

		const double x1 = xmin + step * i;
		const double x2 = xmin + (i + 1) * step;
		#pragma opm atomic 
		simpson_integral += (x2 - x1) / 6.0 * (func(x1) + 4.0 * func(0.5 * (x1 + x2)) + func(x2));
	}

	return simpson_integral;
}


double SimpsonRule::seq_simpson(double (*func)(double x)) {
	double simpson_integral = 0;

	for (long long i = 0; i < n; i ++) {

		const double x1 = xmin + step * i;
		const double x2 = xmin + (i + 1) * step;

		simpson_integral += (x2 - x1) / 6.0 * (func(x1) + 4.0 * func(0.5 * (x1 + x2)) + func(x2));
	}
	
	return simpson_integral;
}