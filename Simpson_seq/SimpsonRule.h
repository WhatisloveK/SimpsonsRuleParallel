#pragma once

 class SimpsonRule {
	 long long n;
	 double xmin, xmax, step;
 public:
	 double seq_simpson(double (*func)(double x));
	 void set_step(double min, double max, long long N);
};
