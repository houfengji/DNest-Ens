#ifndef EXOPLANET_HYPERPARA_H
#define EXOPLANET_HYPERPARA_H

#include <vector>

using namespace std;

struct ExoplanetHyperpara {
	ExoplanetHyperpara();
	vector< vector<double> > uniOrbB;    // uniform prior, Orb for Orbital parameters, B for bounds
	vector<double> uniJtqB;              // uniform prior, B for bounds,  Jitter squared
	vector<double> uniV0B;               // uniform prior, B for bounds
	vector<double> uniV1B;               // uniform prior, B for bounds,  Linear Trend
	
	vector<double> JfOmega;
	vector<double> JfK;
	vector<double> JfJtq;
	
	vector<double> BetaEcc;
};

#endif
