#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "acor.h"
#include "data.h"
#include "evidence.h"
#include "exoplanetjd.h"
#include "int2str.h"
#include "model.h"
#include "preprocessors.h"
#include "quicksort.h"
#include "rng.h"
#include "sampling.h"
#include "sampling_DNest.h"
#include "uniformtest.h"


using namespace std;


void writechain (const Model & model,
                 const vector< vector<double> > & ensemble,
                 bool rep);

int main(void) {
	time_t begin = time(NULL);
	init_genrand(begin);
	time_t end;
	
	//string data_file_name("ploy_3_151_15_-12_9_-1"); size_t dim = 4;
	//string data_file_name("ploy_1_10_5_10"); size_t dim = 2;
	//string data_file_name("uni_2"); size_t dim = 2; size_t ens_size = 20; int num_level = 10; size_t step_size = 10;
	//string data_file_name("uni_10"); size_t dim = 10; size_t ens_size = 50; int num_level = 40; size_t step_size = 100;
	//string data_file_name("122"); size_t num_comp = 1; size_t ens_size = 150; int num_level = 90; size_t step_size = 10;
	string data_file_name("122"); size_t num_comp = 2; size_t ens_size = 200; int num_level = 110; size_t step_size = 10;
	//string data_file_name("122"); size_t num_comp = 3; size_t ens_size = 200; int num_level = 125; size_t step_size = 10; num_level = 1;
	//string data_file_name("gliese_581"); size_t num_comp = 6; size_t ens_size = 200; int num_level = 150; size_t step_size = 10;
	Data data(data_file_name);
	string weight_type = "uniform";
	ExoplanetJD model(data, num_comp, weight_type);
	matrix ensemble;
	model.init(ens_size, ensemble, 0.0000);
	cout << model.LnLikelihood(ensemble[0]) << endl; 
	cout << model.LnDensity(ensemble[0]) << endl;
	
	double ini = 0.00000001;
	size_t burn_in_steps = 5000;
	size_t chain_length  = 50000;
	bool succeed;
	double a = 2;
	
	size_t orbit_crossed = 0;
	
	time_t b = time(NULL);
	model.init(ens_size, ensemble, 0.01);
	for (size_t i = 0; i < burn_in_steps; ++i) {
		sampling(model, ensemble, step_size, a);
	}
	
	for (size_t i = 0; i < chain_length; ++i) {
		sampling(model, ensemble, step_size, a);
		//writechain(model, ensemble, 0);
		
		for (size_t k = 0; k < ens_size; ++k) {
			if (model.orbit_cross(ensemble[k])) {
				++orbit_crossed;
			}
		}
		
	}
	cout << "Orbit Cross: " << orbit_crossed << endl;
	cout << "Total :      " << ens_size*chain_length << endl;
	cout << "Prior Norm Change : " << log (1.0-(double)(orbit_crossed)/(double)(ens_size*chain_length)) << endl;
	
	return 0;
}

void writechain (const Model & model,
                 const vector< vector<double> > & ensemble,
                 bool rep) {

	vector<double> p(model.dim, 0);
	fstream out(("chain_prior_0_" + model.model_name + ".txt").c_str(), ios::out | ios::app);
	double LL;
	for (size_t i = 0; i < ensemble.size(); ++i) {
		//LL = model.LnLikelihood(ensemble[i]);
		if (rep == 1) {
			model.reparametrize(ensemble[i], p);
		}
		else {
			p = ensemble[i];
		}
		for (int64_t j = 0; j < model.dim; ++j) {
			out << setprecision(15) << p[j] << "    ";
		}
		out << endl;
		//if (i != ens_size - 1)
		//	out << endl;
	}
}

