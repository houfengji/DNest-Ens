#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "model.h"
#include "rng.h"

using namespace std;

double sampling (Model & model,      // data and model
                 vector< vector<double> > & ensemble,
                 const size_t num_step,
                 const double a) {         // the tuning in the ensemble sampler

	size_t ens_size = ensemble.size();
	size_t dim      = model.dim;
	
	if (model.dim != ensemble[0].size()) {
		cout << "WARNING: Model dim is not the same as Ensemble dim." << endl;
	}
	
	vector<double> proposed_walker(dim, 0.0);
	size_t choose;
	double random, Z;
	double new_density, old_density;
	double accept;
	size_t accepted = 0;
	
	for (size_t i = 0; i < num_step; ++i) {
		for (size_t k = 0; k < ens_size; ++k) {
			//choose a walker from the complementary ensemble which doesn't include walker_k
			int choose_fail = 0;
			do {
				choose = genrand_int32() % ens_size;
				//if(abs(static_cast<long>(levels[k])-static_cast<double>(levels[choose]))>2&&(++choose_fail)<=100) {
					//continue;
				//}
			} while (choose == k || choose == ens_size);
		
			random = genrand_real2();
			//Z is drawn from a distribution satisfying g(z)=g(1/z)/z.
			//The distribution recommanded in Goodman and Weare's paper is used here.
			//To sample this distribution, direct sampling is the easiest.
			Z = ((a - 1.0) * random + 1.0) * ((a - 1.0) * random + 1.0) / a;
			
			//proposal based on stretch move
			for (size_t j = 0; j < dim; ++j) {
				//X_j(t+1) = Y_j(t) + Z * (X_j(t) - Y_j(t))
				//where Y belongs to the complementary ensemble
				proposed_walker[j] = ensemble[choose][j] + Z * (ensemble[k][j] - ensemble[choose][j]);
				
			}
		
			new_density = model.LnDensity(proposed_walker);
			if (new_density < -1e100) {
				accept = 0;
			}
			else {
				old_density = model.LnDensity(ensemble[k]);
				if (new_density + (dim - 1.0) * log(Z) > old_density) {
					accept = 1;
				}
				else {
					accept = pow(Z, static_cast<int>(dim - 1.0)) * exp(new_density - old_density);
				}
			}
	
			//accept or reject based on accept
			random = genrand_real2();
			if (accept > random) {
				for (size_t j = 0; j < dim; ++j) {
					ensemble[k][j] = proposed_walker[j];
				}
				//ensemble[k] = proposed_walker;
				accepted += 1;
			}
		}
	}
	return (double)accepted/(double)num_step;
}		
