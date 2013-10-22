#include <cmath>
#include <cstdlib>
#include <vector>
#include "exoplanet_hyperpara.h"

using namespace std;

ExoplanetHyperpara::ExoplanetHyperpara(void) {
	uniOrbB.resize(5, vector<double>(2));
	uniOrbB[0][0] = 0;       // lower bound, amp
	uniOrbB[0][1] = 10000;   // upper bound, amp
	uniOrbB[1][0] = 0;       // lower bound, omega
	uniOrbB[1][1] = M_PI;    // upper bound, omega, period 2 day
	uniOrbB[2][0] = 0;       // lower bound, phi
	uniOrbB[2][1] = 2*M_PI;  // upper bound, phi
	uniOrbB[3][0] = 0;       // lower bound, ecc
	uniOrbB[3][1] = 1;       // upper bound, ecc
	uniOrbB[4][0] = 0;       // lower bound, varpi
	uniOrbB[4][1] = 2*M_PI;  // upper bound, varpi
	
	uniJtqB.resize(2);
	uniJtqB[0] = 0;       // lower bound
	uniJtqB[1] = 100000;  // upper bound
	
	uniV0B.resize(2);
	uniV0B[0] = -5000;    // lower bound
	uniV0B[1] = +5000;    // upper bound
	
	JfOmega.resize(1);
	JfOmega[0] = 0.01;
	
	JfK.resize(1);
	JfK[0] = 10;
	
	JfJtq.resize(1);
	JfJtq[0] = 100;
	
	BetaEcc.resize(2);
	BetaEcc[0] = 1;       // A
	BetaEcc[1] = 5;       // B
}
	
