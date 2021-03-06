
#include <iostream>
#include <chrono> 
#include <conio.h>
#include "SimpsonRule.h"
#include<omp.h>
#include <iomanip>

double f(double x) {
	return 2*x + 3*sqrt(x);
}

void output(std::chrono::steady_clock::time_point start, 
	std::chrono::steady_clock::time_point end, double res) {
	double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	std::cout << "result: " << res << " time : " << std::fixed << time_taken << std::setprecision(9) << " sec" << std::endl;
}

int main(int argc, char ** argv)
{
	SimpsonRule simpson_helper;
	simpson_helper.set_step(0, 5000, 3e9);


	auto start = std::chrono::high_resolution_clock::now();
	double res = simpson_helper.seq_simpson(f);
	auto end = std::chrono::high_resolution_clock::now();

	output(start, end, res);

	omp_set_dynamic(false);
	omp_set_num_threads(4);
	start = std::chrono::high_resolution_clock::now();
    res = simpson_helper.parallel_simpson(f);
	end = std::chrono::high_resolution_clock::now();

	output(start, end, res);
	
	return 0;
}

