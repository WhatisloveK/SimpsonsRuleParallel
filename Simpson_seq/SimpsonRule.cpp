#include "SimpsonRule.h"
void SimpsonRule::set_step(double xmin, double xmax, long long n) {
	this->xmax = xmax;
	this->xmin = xmin;
	this->n = n;
	this->step = (xmax - xmin) / n;
}

double simpsonIntegral(double a, double b, int n, const std::function<double(double)>& f) {
	const double width = (b - a) / n;

	double simpson_integral = 0;
	for (int step = 0; step < n; step++) {
		const double x1 = a + step * width;
		const double x2 = a + (step + 1) * width;

		simpson_integral += (x2 - x1) / 6.0 * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
	}

	return simpson_integral;
}
double SimpsonRule::seq_simpson(double (*func)(double x)) {
	double simpson_integral = 0;
	for (long long i = 0; i < n; i ++) {
		//for even numbers
		_sum4 += func(xmin + i * step);
		//for uneven
		_sum2 += func(xmin + (i + 1) * step);
		const double x1 = xmin + step * i;
		const double x2 = xmin + (i + 1) * step;

		simpson_integral += (x2 - x1) / 6.0 * (func(x1) + 4.0 * func(0.5 * (x1 + x2)) + func(x2));
	}
	
	return simpson_integral;
}