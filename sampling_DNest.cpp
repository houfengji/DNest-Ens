#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "evidence.h"
#include "preprocessors.h"
#include "int2str.h"
#include "model.h"
#include "quicksort.h"
#include "randn.h"
#include "rng.h"
#include "sampling_DNest.h"

using namespace std;

double sampling_DNest (Model & model,
                       const size_t num_step,
                       const size_t step_size,
                       vector< vector<double> > & ensemble,
                       vector<size_t> & levels,
                       double a) {
                       
  size_t ens_size = ensemble.size();
	size_t dim      = model.dim;
	
	if (model.dim != ensemble[0].size()) {
		cout << "WARNING: Model dim is not the same as Ensemble dim." << endl;
	}
  
	//model.level.change_weight_type("interval-uniform-e-quantile");
	model.level.change_weight_type("uniform");
	
	
	string out_file_name("post_lev_sample_" + model.model_name + "_" + model.time_label + ".txt");
	fstream out;
	if (num_step == 100000) {
		out.open(out_file_name.c_str(), ios::out | ios::app);
	}

	
	double accept, random;
	double accept_r_i = 0;
	double accept_r_e = 0;
	double accept_r;
	double ens_num_above = 0;
	for (size_t i = 0; i < num_step; ++i) {
		
		accept_r_i = 0;
		accept_r_e = 0;
		for (size_t j = 0; j < step_size; ++j) {
			if (genrand_real2() < 0.5) {
				accept_r_e += update_ensemble(model, ensemble, levels, a);
				accept_r_i += update_index(model, ensemble, levels);
			}
			else {
				accept_r_i += update_index(model, ensemble, levels);
				accept_r_e += update_ensemble(model, ensemble, levels, a);
			}
		}
		//cout << "sampling_DNest - one step: Acceptance Ratio for Index is " << accept_r_i/static_cast<double>(ens_size*step_size) << endl;
		//cout << "sampling_DNest - one step: Acceptance Ratio for Ensem is " << accept_r_e/static_cast<double>(ens_size*step_size) << endl;
		accept_r += 0.5*(accept_r_e + accept_r_i);
		
		if (genrand_real2() < 0.0001) {
			for (size_t k = 0; k < ensemble.size(); ++k) {
				for (size_t l =0 ;l < ensemble[0].size(); ++l) {
					out << setprecision(15) << ensemble[k][l] << "    ";
				}
				out << setprecision(15) << model.LnLikelihood(ensemble[k]) << "     " << levels[k] << endl;
			}
		}
		
		ens_num_above = 0;
		for (size_t k = 0; k < ens_size; ++k) {
			model.chain_LnLikelihood.push_back( model.LnLikelihood(ensemble[k]) );
			++model.total_num_visit;
			++model.Level_Visits_T[levels[k]]; // counting visits
			++model.Level_Visits[levels[k]]; // counting visits
			if (levels[k] < model.level.num_level - 1) {
				if (model.LnLikelihood(ensemble[k]) > model.level.LnThres[levels[k]+1]) {
					++model.Above_Visits[levels[k]]; // counting visits of above
					++ens_num_above;
				}
				if (genrand_real2() < 0.002) {
				cout << "step " << i << " level " << levels[k] << " visits " << model.Level_Visits[levels[k]] << " aboves " << model.Above_Visits[levels[k]] << " acc ratio " << accept_r/static_cast<double>(ens_size*step_size*(i+1)) << endl;
				}
			}
		}
		model.Quant_Ens_Mean.push_back(ens_num_above / static_cast<double>(ens_size));
	}
	
	out.close();
	
	return static_cast<double>(accept_r) / static_cast<double>(num_step);
}
	

// The following routine builds one new level.
// The return value is whether the level is successfully built.
bool new_level (Model & model,
                vector< vector<double> > & ensemble,
                vector< size_t> & levels,
                const double a,
                const size_t burn_in_steps,   // number of steps before actually building the chain
                const size_t chain_length,    // length of the chain used to determine next level
                const size_t chain_length_l,
                const size_t step_size) {     // number of steps between two entries in the chain

	size_t ens_size = ensemble.size();
	size_t dim      = model.dim;
	
	if (model.dim != ensemble[0].size()) {
		cout << "WARNING: Model dim is not the same as Ensemble dim." << endl;
	}
  
  model.level.change_weight_type("exponential");
  
  model.total_num_visit = 0;
  model.Level_Visits_T.clear();
  model.Level_Visits_T.resize(model.level.num_level);
  
	vector<double> LLchain;  // To store the log likelihoods of the chain to find statistics
	
	double tmp;              // Used to determine if its log likelihood is larger than L_star
	
	double accept_r_i = 0;
	double accept_r_e = 0;
	
	for (size_t i = 0; i < burn_in_steps; ++i) {
		for (size_t j = 0; j < step_size; ++j) {
			if (genrand_real2() < 0.5) {
				accept_r_e += update_ensemble(model, ensemble, levels, a);
				accept_r_i += update_index(model, ensemble, levels);
			}
			else {
				accept_r_i += update_index(model, ensemble, levels);
				accept_r_e += update_ensemble(model, ensemble, levels, a);
			}
		}
		for (size_t k = 0; k < ens_size; ++k) {
			++model.total_num_visit;
			++model.Level_Visits_T[levels[k]]; // counting visits
			if ((double)rand() < 0.0001*(double)RAND_MAX) {
				cout << "step " << i << " out of " << burn_in_steps << " level " << levels[k] << " visits " << model.Level_Visits_T[levels[k]] << endl;
			}
		}
	}
	
	//cout << "build_level - burning: Acceptance Ratio for Index is " << accept_r_i/(step_size*burn_in_steps*ens_size) << endl;
	//cout << "build_level - burning: Acceptance Ratio for Ensem is " << accept_r_e/(step_size*burn_in_steps*ens_size) << endl;
	
	bool chain_length_done = 0;
	while (LLchain.size()<chain_length_l) {
		size_t old_length = LLchain.size();
		accept_r_i = 0;
		accept_r_e = 0;
		for (size_t j = 0; j < step_size; ++j) {
			if (genrand_real2() < 0.5) {
				accept_r_e += update_ensemble(model, ensemble, levels, a);
				accept_r_i += update_index(model, ensemble, levels);
			}
			else {
				accept_r_i += update_index(model, ensemble, levels);
				accept_r_e += update_ensemble(model, ensemble, levels, a);
			}
		}
		//cout << "build_level - build chain one step: Acceptance Ratio for Index is " << accept_r_i/(step_size*ens_size) << endl;
		//cout << "build_level - build chain one step: Acceptance Ratio for Ensem is " << accept_r_e/(step_size*ens_size) << endl;
		
		
		
		for (size_t k = 0; k < ens_size; ++k) {
			
			tmp = model.LnLikelihood(ensemble[k]);
			
			if (tmp > model.level.LnThres.back()) {
				LLchain.push_back(tmp);
			}
		}
		
		size_t new_length = LLchain.size();
		//if (new_length - old_length < ens_size/10) {
			//cout << new_length << endl;
		//}
		if (new_length == old_length) {
			cerr << "No Samples Above Latest Level Found!" << endl;
			return 0;
		}
		// Randomly monitor the code.
		//if (genrand_real2() < 0.1) {
			//cout << "Length of Chain = " << new_length << endl;
		//}
		for (size_t k = 0; k < ens_size; ++k) {
			++model.total_num_visit;
			++model.Level_Visits_T[levels[k]]; // counting visits
			if ((double)rand() < 0.0001*(double)RAND_MAX) {
				cout << "Length of Chain = " << new_length << " out of " << chain_length << " level " << levels[k] << " visits " << model.Level_Visits_T[levels[k]] << endl;
			}
		}
		
		if (LLchain.size() > chain_length && chain_length_done == 0) {
			chain_length_done = 1;
			for (size_t i = 0; i < model.Level_Visits_T.size(); ++i) {
				cout << "Level " << i << " Visits " << model.Level_Visits_T[i] << endl;
			}
			LLchain.resize(chain_length-1);
	
			quicksort(LLchain, LLchain.begin(), LLchain.end()-1);  // rank likelihoods in descending order
			size_t new_lnthres_index = static_cast<size_t>((1.0+LLchain.size())/exp(1)) - 1;
			if ( fabs((double)(LLchain.size()+1)/(double)(new_lnthres_index+1)-exp(1)) > fabs((double)(LLchain.size()+1)/(double)(new_lnthres_index+2)-exp(1)) ) {
				new_lnthres_index = new_lnthres_index + 1;
			}
			//cout << new_llevel_index << endl;
			vector<double> inbetween(LLchain.size() - new_lnthres_index);
			for (size_t k = 0; k < inbetween.size(); ++k) {
				inbetween[k] =  LLchain[LLchain.size()-1-k];
			}
			double new_lnthres = LLchain[new_lnthres_index];
			if (new_lnthres < -1e100) {
				for (size_t i = 0; i < LLchain.size(); ++i) {
					if (LLchain[i] < -1e100) {
						new_lnthres_index = i-1;
						new_lnthres = LLchain[new_lnthres_index];
						break;
					}
				}
			}
			double new_lnprims = model.level.LnPrims.back() + log(new_lnthres_index + 1) - log(LLchain.size() + 1);
			size_t n = model.level.num_level;
	
			if (n>=4 && new_lnthres - model.level.LnThres[n-1] > model.level.LnThres[n-1] - model.level.LnThres[n-2]) {
				;
			}
			else {
				model.level.add_level(new_lnthres, new_lnprims);
				return 1;
			}
		}
	}
	
	
	for (size_t i = 0; i < model.Level_Visits_T.size(); ++i) {
		cout << "Level " << i << " Visits " << model.Level_Visits_T[i] << endl;
	}
	LLchain.resize(chain_length-1);
	
	quicksort(LLchain, LLchain.begin(), LLchain.end()-1);  // rank likelihoods in descending order
	size_t new_lnthres_index = static_cast<size_t>((1.0+LLchain.size())/exp(1)) - 1;
	if ( fabs((double)(LLchain.size()+1)/(double)(new_lnthres_index+1)-exp(1)) > fabs((double)(LLchain.size()+1)/(double)(new_lnthres_index+2)-exp(1)) ) {
		new_lnthres_index = new_lnthres_index + 1;
	}
	
	//cout << new_llevel_index << endl;
	vector<double> inbetween(LLchain.size() - new_lnthres_index);
	for (size_t k = 0; k < inbetween.size(); ++k) {
		inbetween[k] =  LLchain[LLchain.size()-1-k];
	}
	double new_lnthres = LLchain[new_lnthres_index];
	if (new_lnthres < -1e100) {
		for (size_t i = 0; i < LLchain.size(); ++i) {
			if (LLchain[i] < -1e100) {
				new_lnthres_index = i-1;
				new_lnthres = LLchain[new_lnthres_index];
				break;
			}
		}
	}
	
	double new_lnprims = model.level.LnPrims.back() + log(new_lnthres_index + 1) - log(LLchain.size() + 1);
	model.level.add_level(new_lnthres, new_lnprims);
	
	return 1;
}

double LLmean(Model & model, matrix & ensemble) {
	size_t ens_size = ensemble.size();
	double LL = 0;
	for (size_t i = 0; i < ens_size; ++i) {
		LL += model.LnLikelihood(ensemble[i]);
	}
	return LL/static_cast<double>(ens_size);
}

void ensemble_mean(const matrix & ensemble, std::vector<double> & mean) {
	size_t ens_size = ensemble.size();
	size_t dim = ensemble[0].size();
	if (dim != mean.size()) {
		mean.resize(dim);
	}
	for (size_t i = 0; i < dim; ++i) {
		mean[i] = 0;
		for (size_t j = 0; j < ens_size; ++j) {
			mean[i] += ensemble[j][i];
		}
		mean[i] /= static_cast<double>(ens_size);
	}
}

double accept_ensemble (Model & model, 
                        const size_t & level,
                        const matrix & ensemble, 
                        const matrix & proposed_ensemble) {
	size_t ens_size = ensemble.size();
	size_t dim = ensemble[0].size();
	double new_density = 0;
	double old_density = 0;
	for (size_t k = 0; k < ens_size; ++k) {
		new_density += model.LnDensity( proposed_ensemble[k], level);
		old_density += model.LnDensity( ensemble[k], level);
	}
	new_density /= static_cast<double>(ens_size);
	old_density /= static_cast<double>(ens_size);
	
	double accept = exp(new_density) * exp(-old_density);
	return accept;
}

// With different walkers belonging to different levels,
// the following routine updates the whole ensemble one step.
// Returned value is the number of moves accepted in the routine.
size_t update_ensemble (Model & model,      // data and model
                        vector< vector<double> > & ensemble,
                        const vector<size_t> & levels,
                        const double a) {         // the tuning in the ensemble sampler

	size_t ens_size = ensemble.size();
	size_t dim      = model.dim;
	
	if (model.dim != ensemble[0].size()) {
		cout << "WARNING: Model dim is not the same as Ensemble dim." << endl;
	}
	if (ens_size != levels.size()) {
		cerr << "ERROR: Number of levels is not the same as Ensemble size!" << endl;
		return 0;
	}
	
	vector<double> proposed_walker(dim, 0.0);
	size_t choose;
	double random, Z;
	double new_density, old_density;
	double accept;
	size_t accepted = 0;
	
	for (size_t k = 0; k < ens_size; ++k) {
		//choose a walker from the complementary ensemble which doesn't include walker_k
		int choose_fail = 0;
		do {
			choose = genrand_int32() % ens_size;
			//if(abs(static_cast<long>(levels[k])-static_cast<long>(levels[choose]))>2&&(++choose_fail)<=100) {
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
		
		new_density = model.LnDensity(proposed_walker, levels[k]);
		if (new_density < -1e100) {
			accept = 0;
		}
		else {
			old_density = model.LnDensity(ensemble[k], levels[k]);
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
	return accepted;
}

// With different walkers belonging to different levels,
// the following routine updates the level index one step.
// Returned value is the number of steps accepted in the routine.
size_t update_index (Model & model,
                     const vector< vector<double> > & ensemble,
                     vector<size_t> & levels) {
                     
	size_t ens_size = ensemble.size();
	size_t dim      = model.dim;
	
	if (model.dim != ensemble[0].size()) {
		cout << "WARNING: Model dim is not the same as Ensemble dim." << endl;
	}
	if (ens_size != levels.size()) {
		cerr << "ERROR: Number of levels is not the same as Ensemble size!" << endl;
		return 0;
	}
	
	double random;
	double new_density, old_density, new_weight, old_weight;
	double accept;
	size_t accepted = 0;
  size_t proposed_level;
  long _proposed_level;
  
  for (size_t k = 0; k < ens_size; ++k) {
		
		if (levels[k] > 0 && levels[k] < model.level.num_level-1) {
			proposed_level = (genrand_real2() < 0.5) ? (levels[k] + 1) : (levels[k] - 1);
		}
		else {
			if (levels[k] == 0) {
				proposed_level = (genrand_real2() < 0.5) ? (levels[k] + 1) : (levels[k]);
				if (model.level.num_level == 1) {
					proposed_level = levels[k];
				}
			}
			else {
				proposed_level = (genrand_real2() < 0.5) ? (levels[k] - 1) : (levels[k]);
			}
		}
		
		new_weight = model.level.LWeight[proposed_level];
		old_weight = model.level.LWeight[levels[k]];
		
		if (proposed_level > levels[k]) {
	  	if (model.LnDensity(ensemble[k], proposed_level) < -1e100 ) {
	  		accept = 0;
	  	}
	  	else {
	  		accept = 1;
	  	}
		}
		else {
			new_density = model.LnDensity(ensemble[k], proposed_level);
			old_density = model.LnDensity(ensemble[k], levels[k]);
			
			accept = exp(new_density - old_density) * exp(new_weight - old_weight);
		}
		
		// THE FOLLOWING IF BLOCK VIOLATES THE MARKOV PROPERTY!
		
		if (model.Level_Visits_T.size() == model.level.num_level) {
			double n_ori = model.Level_Visits_T[levels[k]];        // the actual number of visits to the original level
			double n_pro = model.Level_Visits_T[proposed_level];   // the acutal number of visits to the proposed level
			double n_ori_theo = model.total_num_visit * exp(model.level.LWeight[levels[k]]-model.level.LWeight_norm);  // the theoretical number of visits to the original level
			double n_pro_theo = model.total_num_visit * exp(model.level.LWeight[proposed_level]-model.level.LWeight_norm);  // the theoretical number of visits to the proposed level
			double C = 10.0;
			double beta = 500.0;
			double alpha = pow((n_ori+C)*(n_pro_theo+C)/(n_ori_theo+C)/(n_pro+C),beta);
			if (isinf(alpha)) {
				accept = 1;
			}
			else {
				accept *= alpha;  // correction to the acceptance probability
			}
		}
		
		
		random = genrand_real2();
		if (accept > random) {
			levels[k] = proposed_level;
			accepted += 1;
		}
	}
	
	return accepted;
}
