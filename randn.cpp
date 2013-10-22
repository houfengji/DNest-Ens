/*
 *  This cpp file defines the function to generate 1-D standard
 *  normal i.i.d. using Box Muller.
 */
 
#include <cmath>
#include "rng.h"
#include "randn.h"

double randn ( void ) {
	double rand1 = genrand_real2();
	double rand2 = genrand_real2();
	double R = sqrt(-2.0 * log(rand1));
	double T = 2 * M_PI * rand2;
	return R * cos(T);
}
