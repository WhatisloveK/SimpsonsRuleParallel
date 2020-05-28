
#include <iostream>
#include <ctime>
#include <conio.h>
#include "SimpsonRule.h"
#include<omp.h>

double f(double x) {
	return 2*x + 3*sqrt(x);
}

int main(int argc, char ** argv)
{
	SimpsonRule simpson_helper;
	simpson_helper.set_step(0, 5, 2e8);
	time_t start, finish;
	double duration;

	/*srand(time(0));
	int count = 1000000;

	int *ptrarray = new int[count];
	for (int i = 0; i < count; i++)
		ptrarray[i] = rand() % 10000;*/

	start = clock();
	
	double res = simpson_helper.seq_simpson(f);

	finish = clock();
	duration = (finish - start) / double(CLOCKS_PER_SEC);
	std::cout << "result: " << res << " time: " << duration<<std::endl;


	omp_set_num_threads(4);
	double tm = omp_get_wtime();
	res = simpson_helper.parallel_simpson(f);
	tm = omp_get_wtime() - tm;

	
	
	std::cout <<"result: "<<res<<" time: "<< tm;

	system("pause");
	return 0;
}

