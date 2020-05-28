#include "SimpsonRule.h"
void SimpsonRule::set_step(double xmin, double xmax, long long n) {
	this->xmax = xmax;
	this->xmin = xmin;
	this->n = n;
	this->step = (xmax - xmin) / n;
}

double SimpsonRule::seq_simpson(double (*func)(double x)) {
	double _sum4 = 0, _sum2 = 0;
	for (long long i = 1; i < 2 * n; i += 2) {
		//for even numbers
		_sum4 += func(xmin + i * step);
		//for uneven
		_sum2 += func(xmin + (i + 1) * step);
	}
	double res = step / 3 * (func(xmin) + 4 * _sum4 + 2 * _sum2 - func(xmax));
	return res;
}