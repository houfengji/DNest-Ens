/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *
 */
 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "data.h"
#include "preprocessors.h"
#include "int2str.h"
#include "level.h"
#include "model.h"
#include "rng.h"
#include "uniformtest.h"

using namespace std;


// Constructor
UniformTest::UniformTest(Data & dt,                 // Data for likelihood calculation
                         size_t dimension,          // Dimension of the model
                         std::string weight_type):  // Name for the weight type assigned to each level
                         data(dt),
                         Model(dimension, weight_type) {
	model_name = "uniformtest_" + int2str(dim,3);
	upper_bound = 10;
	lower_bound = -10;
	best_fit.resize(dim+1);
	best_fit[dim] = -1e300;
}

// Reparametrization (not used so far)
double UniformTest::reparametrize(const vector<double> & param, vector<double> & p) const
{
	for (int64_t i = 0; i < dim; ++i) {
		p[i] = param[i];
	}
	return 1;
}

// Log Prior
// All the parameters have Uniform Prior
double UniformTest::LnPrior(const vector<double> & param) 
{
	vector<double> p(dim, 0);
	reparametrize(param, p);
	
	double lp = 0;
	for (size_t i = 0; i < dim; ++i) {
		if (p[i] < lower_bound || p[i] > upper_bound) {
			return -1e300;
		}
		lp += - log(upper_bound-lower_bound);
	}
	
	return lp;
}

// Log Prior with Threshold
double UniformTest::LnPrior(const vector<double> & param,
                            const size_t j)           // the index of the level
{
	if (j >= level.num_level) {
		cerr << "UniformTest, LnPrior: level index out of range!" << endl;
	}
	double LL_star = level.LnThres[j];  // current log likelihood threshold
	double LX_star = level.LnPrims[j];  // corresponding X-Value
	double LL = LnLikelihood(param);   // log likelihood of the current parameters
	if (LL < LL_star) {
		return -1e308;      // When the current log likelihood is smaller than the threshold,
		                    // return an extremely small number.
	}
	else {
		vector<double> p(dim, 0);
		reparametrize(param, p);
		
		/*
		for (size_t i = 0; i < dim-1; ++i) {
			if (p[i] < p[i+1] ) {
				return -1e300;
			}
		}
	*/
	
		double lp = 0;
		for (size_t i = 0; i < dim; ++i) {
			if (p[i] < lower_bound || p[i] > upper_bound) {
				return -1e300;
			}
			lp += - log(upper_bound-lower_bound);
		}
		return lp - LX_star;   // X-Value is used a normalization term
	}
}

// Log Likelihood
// Gaussian Likelihood
double UniformTest::LnLikelihood(const vector<double> & param)
{
	vector<double> p(dim, 0);
	reparametrize(param, p);
	
	double ln_error_sq = 0;
	double chi_square = 0;
	double poly_value = 0;
	for (int64_t i = 0; i < data.data_size; ++i) {
	
		for (int64_t j = 0; j < dim; ++j) {	
			chi_square += pow((data.data[i][j] - p[j]), 2) / pow(data.data[i][dim], 2);
		}
		ln_error_sq += 2 * log(data.data[i][dim]);
	}
	return - 0.5 * static_cast<double>(dim) * data.data_size * log(2 * M_PI) - 0.5 * static_cast<double>(dim) * ln_error_sq - 0.5 * chi_square;
}
 
double UniformTest::LnDensity(const vector<double> & param) 
{
	return LnPrior(param);
}

double UniformTest::LnDensity(const vector<double> & param, const size_t j) 
{
	return LnPrior(param, j);
}

// initialize the ensemble
void UniformTest::init(size_t ens_size, matrix & ensemble, double ini) const
{
	ensemble.resize(ens_size, vector<double>(dim, 0));
	if (data.data_name.compare("uni_1") == 0 && dim == 1) {
		for (size_t i = 0; i < ens_size; ++i) {
			ensemble[i][0] = 0.01 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
		}
	}
	if (data.data_name.compare("uni_2") == 0 && dim == 2) {
		for (size_t i = 0; i < ens_size; ++i) {
			ensemble[i][0] = 0.01 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][1] = -0.01 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
		}
	}
	if (data.data_name.compare("uni_10") == 0 && dim == 10) {
		for (size_t i = 0; i < ens_size; ++i) {
			ensemble[i][0] = 0.05 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][1] = 0.04 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][2] = 0.03 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][3] = 0.02 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][4] = 0.01 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][5] = -0.01 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][6] = -0.02 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][7] = -0.03 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][8] = -0.04 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			ensemble[i][9] = -0.05 * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
		}
	}
	if (data.data_name.compare("uni_100") == 0 && dim == 100) {
		for (size_t i = 0; i < ens_size; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				ensemble[i][j] = 0.000001 * (((double)rand() / (double)RAND_MAX - 0.5) * 2) * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			}
		}
	}
}

UniformTest::~UniformTest() {
}
	
	
double UniformTest::theo_lnprims_2d(double ln_likelihood) {
	double area = - 2*M_PI*log(2*M_PI) - 2*M_PI*ln_likelihood;
	return log(area/400.);
}
