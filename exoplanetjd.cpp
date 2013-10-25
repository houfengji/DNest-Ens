/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  This cpp file defines all the member functions in ExoplanetJD class.
 */
 
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "data.h"
#include "int2str.h"
#include "keplers_eqn.h"
#include "level.h"
#include "model.h"
#include "exception.h"
#include "exoplanet_hyperpara.h"
#include "exoplanet_init.h"
#include "exoplanetjd.h"

using namespace std;


size_t ExoplanetJD::count_distinction(void) {
	if (data.num_col == 3) {
		return 1; // If there are 3 columns, then the data is single-sourced.
	}
	if (data.num_col == 4) { // The 4th column is usually used for distinctions
		vector<size_t> dist(data.num_row, 0);
		for (size_t i = 0; i < data.num_row; ++i) {
			dist[i] = static_cast<size_t>(data.data[i][3]);
		}
		
		vector<size_t> distinctions;
		distinctions.push_back(dist[0]);
		bool found;
		for (size_t i = 1; i < data.num_row; ++i) {
			found = 0;
			for (size_t j = 0; j < distinctions.size(); ++j) {
				if (dist[i] == distinctions[j]) {
					found = 1;
					break;
				}
			}
			if (found == 1) {
				continue;
			}
			else {
				distinctions.push_back(dist[i]);
			}
		}
		
		for (size_t i = 0; i < data.num_row; ++i) {
			found = 0;
			for (size_t j = 0; j < distinctions.size(); ++j) {
				if (dist[i] == distinctions[j]) {
					found = 1;
					data.data[i][3] = j; // re-label the distinctions as 0, 1, 2...
					break;
				}
			}
		}
		return distinctions.size();
	}
	
	cerr << "To use ExoplanetJD, the data is at least required to have 3 or 4 columns." << endl;
	cerr << "The current data set Doesn't satisfy this least requirement." << endl;
	return 0;
}

// Constructor
ExoplanetJD::ExoplanetJD(Data & dt,                 // Data for likelihood calculation
                         size_t num_companion,      // Dimension of the model
                         std::string weight_type):  // Name for the weight type assigned to each level
                         data(dt),
                         num_comp(num_companion),
                         Model(weight_type) {
	double begin_time = data.data[0][0];
	for (size_t i = 0; i < data.data_size; ++i) {
		data.data[i][0] -= begin_time;
	}
	model_name = "exop_" + data.data_name + "_mod_" + int2str(num_comp);
	num_d = count_distinction();
	if (num_d == 0) {
		throw( Exception("Data format can not be processed!") );
	}
	if (num_d > 1) {
		model_name = model_name + "_d";
	}
	dim = 5 * num_comp + 2 * num_d;
	best_fit.resize(dim+1);
	best_fit[dim] = -1e300;
	LnPermutation = 0;
	for (size_t i = 1; i < num_comp; ++i) {
		LnPermutation += log(static_cast<double>(num_comp));
	}
}

// Reparametrization
// Good re-parametrization may improve sampling greatly. But in the context of
// DNS, re-parametrization is not very necessary.
double ExoplanetJD::reparametrize(const vector<double> & param, vector<double> & p) const
{
	//This set of reparametrization was used when we used MCMC for inference.
	//param_0 = omega
	//param_1 = sqrt(amplitude) * sin(phi)
	//param_2 = sqrt(amplitude) * cos(phi)
	//param_3 = sqrt(eccentricity) * sin(varpi)
	//param_4 = sqrt(eccentricity) * cos(varpi)
	
	//This reparametrization has Jacobian 1.
	
	for (size_t i = 0; i < dim; ++i) {
		p[i] = param[i];
	}
	return 1;
}

// This routine checks if the orbits of two planets cross.
// If they do cross, return 1.
// If they do NOT cross, return 0.
bool ExoplanetJD::orbit_cross (const vector<double> & parameter) 
{  // parameters used in MCMC
	if (num_comp <= 1) {
		return 0;
	}

	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);

	double P_large, P_small;
	double a_large, a_small;
	double e_large, e_small;
	for (int64_t i = 0; i < num_comp; ++i) {
		for (int64_t j = i + 1; j < num_comp; ++j) {
			if (p[i*5+1] < p[j*5+1]) {
				P_large = 2 * M_PI / p[i*5+1];
				P_small = 2 * M_PI / p[j*5+1];
				a_large = pow(P_large*P_large, 1.0/3.0);
				a_small = pow(P_small*P_small, 1.0/3.0);
				e_large = p[i*5+3];
				e_small = p[j*5+3];
				if (a_small * (1 + e_small) > a_large * (1 - e_large)) {
					return 1;
				}
			}
			else {
				P_large = 2 * M_PI / p[j*5+1];
				P_small = 2 * M_PI / p[i*5+1];
				a_large = pow(P_large*P_large, 1.0/3.0);
				a_small = pow(P_small*P_small, 1.0/3.0);
				e_large = p[j*5+3];
				e_small = p[i*5+3];
				if (a_small * (1 + e_small) > a_large * (1 - e_large)) {
					return 1;
				}
			}
		}
	}
	return 0;
}

// We force the periods of the companions to be in decreasing order (omega in  ascending order).
// If the periods of the companions are not in ascending order, return 1.
bool ExoplanetJD::exchange (const vector<double> & parameter) 
{
	if (num_comp <= 1) {
		return 0;
	}

	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);

	for (size_t i = 0; i < num_comp-1; ++i) {
		if (p[i*5+1] < p[(i+1)*5+1]) {
			return 1;
		}
	}
	return 0;
}

// This function serves as a sanity check, because we may need to make 
// likelihood calls when the parameters are not physical, i.e. eccentricity > 1.
// If we don't make such check, likelihood may behave humorously. 
// If combinations of parameters are not physical, return 0.
bool ExoplanetJD::physical (const std::vector<double> & parameter)
{
	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	bool crossed = orbit_cross(parameter);
	//forbid orbit crossing
	if (crossed == 1) {
		return 0;  // 0 for unphysical
	}
	for (size_t i = 0; i < num_comp * 5; i=i+5) {
		if (p[i+3] < 0 || p[i+3] >=1) {
			return 0;
		}
		if (p[i+2] < 0) {
			return 0;
		}
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[5*num_comp+num_d+i] < 0) {
			return 0;
		}
	}
	return 1;
}

// All uniform priors
// Which prior is to be used is determined in LnPrior.
// NOTE: Orbit crossing and companion permutation are not included.
double ExoplanetJD::LnPrior01 (const vector<double> & parameter)  // All Uniform Dummy Prior 
{

	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	double lpr = 0;
	for (size_t i = 0; i < num_comp; ++i) {
		for (size_t j = 0; j < 5; ++j) {  // hard-wired FIVE because there are 5 orbital parameters.
			if (p[i*5+j] < hyper.uniOrbB[j][0] || p[i*5+j] > hyper.uniOrbB[j][1]) {
				return -1e300;
			}
			else {
				lpr += - log(hyper.uniOrbB[j][1] - hyper.uniOrbB[j][0]);
			}
		}
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+i] < hyper.uniV0B[0] || p[num_comp*5+i] > hyper.uniV0B[1]) {
			return -1e300;
		}
		else {
			lpr += - log(hyper.uniV0B[1] - hyper.uniV0B[0]);
		}
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+num_d+i] < hyper.uniJtqB[0] || p[num_comp*5+num_d+i] > hyper.uniJtqB[1]) {
			return -1e300;
		}
		else {
			lpr += - log(hyper.uniJtqB[1] - hyper.uniJtqB[0]);
		}
	}
	return lpr;
}

// More sensible and informative priors
// Which prior is to be used is determined in LnPrior.
// NOTE: Orbit crossing and companion permutation are not included.
double ExoplanetJD::LnPrior02 (const vector<double> & parameter)  
{

	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	double lpr = 0;
	for (size_t i = 0; i < num_comp; ++i) {
		for (size_t j = 0; j < 5; ++j) {  // hard-wired FIVE because there are 5 orbital parameters.
			if (p[i*5+j] < hyper.uniOrbB[j][0] || p[i*5+j] > hyper.uniOrbB[j][1]) {
				return -1e300;
			}
		}
		lpr += - log(p[i*5+1] + hyper.JfOmega[0]);
		lpr += - log(p[i*5+0] + hyper.JfK[0]);
		lpr += (hyper.BetaEcc[0]-1)*log(p[i*5+3]) + (hyper.BetaEcc[1]-1)*log(1-p[i*5+3]); 
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+i] < hyper.uniV0B[0] || p[num_comp*5+i] > hyper.uniV0B[1]) {
			return -1e300;
		}
	}
	for (size_t i = 0; i < num_d; ++i) {
		if (p[num_comp*5+num_d+i] < hyper.uniJtqB[0] || p[num_comp*5+num_d+i] > hyper.uniJtqB[1]) {
			return -1e300;
		}
		lpr += - log(p[num_comp*5+num_d+i] + hyper.JfJtq[0]);
	}
	return lpr;
}

// Log Prior with Threshold
// Determines which prior to be used
// Orbit crossing and companion permutation are included.
double ExoplanetJD::LnPrior(const vector<double> & parameter,  // All Uniform Dummy Prior 
                            const size_t j)          // the index of the level
{
	if (j >= level.num_level) {
		throw( Exception("ExoplanetJD, LnPrior: level index out of range!") );
	}
	if ( exchange(parameter) ) {
		return -1e300;
	}
	if ( !physical(parameter) ) {
		return -1e300;   
	}
	double LL_star = level.LnThres[j];  // current log likelihood threshold
	double LX_star = level.LnPrims[j];  // corresponding X-Value
	double LL = LnLikelihood(parameter);   // log likelihood of the current parameters
	if (LL < LL_star) {
		return -1e300;      // When the current log likelihood is smaller than the threshold,
		                    // return an extremely small number.
	}
	else {
		vector<double> p(dim, 0);
		double jacobian = reparametrize(parameter, p);
	
		double lpr = LnPrior02(parameter);      // NOTE: choose prior here
		return lpr - LX_star;   // X-Value is used a normalization term
	}
}
	
	
// the logarithm of likelihood function
double ExoplanetJD::LnLikelihood(const vector<double> & parameter) 
{
	vector<double> p(dim, 0);
	double jacobian = reparametrize(parameter, p);
	// If the eccentricity steps into absurd zones,
  // directly return a small value that means impossible.
	if (physical(parameter) == 0) {  // unphysical
		return -1.e290; // small
	}

	double chi_sq      = 0.0;
	double ln_error_sq = 0.0;
	double rv_pred     = 0.0;
	double ln_jacobian = 0.0;
	double jitter_sq;
	
	for (size_t i = 0; i < data.data_size; ++i) {
		for (size_t j = 0; j < num_comp * 5; j=j+5) {
			rv_pred += rad_v_pred(data.data[i][0], p[j+0], p[j+1], p[j+2], p[j+3], p[j+4]);
			if (isnan(rv_pred) || isinf(rv_pred)) {
				cerr << "Keplers Equation returns unphysical results not understood." << endl;
				cerr << setprecision(15) << data.data[i][0] << "  " << p[j+0] << "  " << p[j+1] << "  " << p[j+2] << "  " << p[j+3] << "  " << p[j+4] << endl;
				cerr << "isnan = " << isnan(rv_pred) << ", isinf = " << isinf(rv_pred) << endl;
				return -1.e290;
			}
		}
		if (data.num_col == 4) {
			rv_pred += p[num_comp*5+static_cast<size_t>(data.data[i][3])];
			jitter_sq = p[num_comp*5+num_d+static_cast<size_t>(data.data[i][3])];
		}
		else {
			if (data.num_col == 3) {
				rv_pred += p[num_comp*5];
				jitter_sq = p[num_comp*5+num_d];
			}
			else {
				cerr << "The data is required to have 3 or 4 columns." << endl;
				cerr << "Zero Likelihood will be returned." << endl;
				return -1e290;
			}
		}
		chi_sq += (data.data[i][1] - rv_pred) * (data.data[i][1] - rv_pred) / (data.data[i][2] * data.data[i][2] + jitter_sq);  // the exponent
    
		ln_error_sq += log(2.0*M_PI*(data.data[i][2] * data.data[i][2] + jitter_sq)); // the normalization
		rv_pred = 0;
	}

	ln_jacobian = log(jacobian);
	
	
	double lnlikelihood = - 0.5 * ln_error_sq - 0.5 * chi_sq - ln_jacobian;
	if (lnlikelihood > best_fit.back()) {
		for (size_t k = 0; k < dim; ++k) {
			best_fit[k] = parameter[k];
		}
		best_fit[dim] = lnlikelihood;
		//cout << setprecision(15) << best_fit[dim] << endl;
		fstream outbest( ("best_fit_"+model_name+"_"+time_label+".txt").c_str(), ios::out);
		for (size_t k = 0; k < dim+1; ++k) {
			outbest << "init[" << k << "] = ";
			outbest << setprecision(15) << best_fit[k] << ";" <<endl;
		}
		outbest.close();
	}
	return - 0.5 * ln_error_sq - 0.5 * chi_sq - ln_jacobian;
}

double ExoplanetJD::LnDensity(const vector<double> & parameter)
{
	if (exchange(parameter) == 1) {
		return -1e300;
	}
	double lpr = LnPrior02(parameter);
	return lpr;   // X-Value is used a normalization term
}

double ExoplanetJD::LnDensity(const vector<double> & parameter, const size_t j)
{
	return LnPrior(parameter, j);
}

void ExoplanetJD::init(size_t ens_size, std::vector< std::vector<double> > & ensemble, double ini)
{
	bool fail = 0;
	srand(time(NULL));
	//122 or any other no dis
	if (num_d == 1) {
		ensemble.resize(ens_size, vector<double>(dim, 0));
		ExoplanetInit initialize(data.data_name, num_comp);
		for (size_t i = 0; i < ens_size; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				ensemble[i][j] = initialize.init[j] * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			}
	
			if (LnLikelihood(ensemble[i]) < level.LnThres.back()) {
				fail = 1;
				break;
			}
		}
		
	}
	
	//282 no dis
	if (data.data_name.compare("282") == 0 && num_d == 1) {
		ensemble.resize(ens_size, vector<double>(dim, 0));
		ExoplanetInit initialize(data.data_name, num_comp);
		for (size_t i = 0; i < ens_size; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				ensemble[i][j] = initialize.init[j] * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			}
	
			if (LnLikelihood(ensemble[i]) < level.LnThres.back()) {
				fail = 1;
				break;
			}
		}
		
	}
	
	//gliese 581 2 dis
	if (data.data_name.compare("gliese_581") == 0 && num_d == 2) {
		ensemble.resize(ens_size, vector<double>(dim, 0));
		ExoplanetInit initialize(data.data_name, num_comp, num_d);
		for (size_t i = 0; i < ens_size; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				ensemble[i][j] = initialize.init[j] * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);
			}
	
			if (LnLikelihood(ensemble[i]) < level.LnThres.back()) {
				fail = 1;
				break;
			}
		}
		
	}
	
	
	if (fail == 1) {
		for (size_t i = 0; i < ens_size; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				ensemble[i][j] = best_fit[j] * (1.0 + ini * ((double)rand() / (double)RAND_MAX - 0.5) * 2);;
			}
		}
	}
}

ExoplanetJD::~ExoplanetJD() {
}
