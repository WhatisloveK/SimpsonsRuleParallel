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

	double sum2 = 0, sum4 = 0, sum = 0;
	double _xmin = xmin, _step = step;
	long long i;
#pragma omp parallel for reduction(+: sum4, sum2) private(i) shared(_xmin,_step)
	for (i = 1; i <= n; i += 2)
	{
		sum4 += func(_xmin + _step * i);//Значения с нечётными индексами, которые нужно умножить на 4.
		sum2 += func(_xmin + _step * (i + 1));//Значения с чётными индексами, которые нужно умножить на 2.
	}
	sum = func(xmin) + 4 * sum4 + 2 * sum2 - func(xmax);//Отнимаем значение f(b) так как ранее прибавили его дважды. 
	simpson_integral = (step / 3) * sum;

	return simpson_integral;
}


double SimpsonRule::seq_simpson(double (*func)(double x)) {
	double simpson_integral = 0;

	double sum2 = 0, sum4 = 0, sum = 0;
	
	for (long long i = 1; i <= n; i += 2)
	{
		sum4 += func(xmin + step * i);//Значения с нечётными индексами, которые нужно умножить на 4.
		sum2 += func(xmin + step * (i + 1));//Значения с чётными индексами, которые нужно умножить на 2.
	}
	sum = func(xmin) + 4 * sum4 + 2 * sum2 - func(xmax);//Отнимаем значение f(b) так как ранее прибавили его дважды. 
	simpson_integral = (step / 3) * sum;
	
	return simpson_integral;
}