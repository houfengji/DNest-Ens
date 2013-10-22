/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  Aug 5, 2013
 *  This cpp file defines a simple constructor for Model Base class.
 *
 */
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include "int2str.h"
#include "model.h"


// Constructor:
// the dimension of the model's paramter space
// weight type for the level
Model::Model(size_t dimension, std::string weight_type): 
             dim(dimension), level(weight_type){

	time_label = int2str(static_cast<long>(time(NULL)));
}

// Constructor:
// weight type for the level
Model::Model(std::string weight_type): 
             level(weight_type){
	// dimension will be given inside the derivative constructor
	time_label = int2str(static_cast<long>(time(NULL)));
}

void Model::clear_all_visits(void) {
	Level_Visits_T.clear();
	Level_Visits.clear();
	Above_Visits.clear();
	Quant_Ens_Mean.clear();
	chain_LnLikelihood.clear();
	total_num_visit = 0;
}

void Model::clear_for_sampling(void) {
	Level_Visits.clear();
	Above_Visits.clear();
	Quant_Ens_Mean.clear();
	chain_LnLikelihood.clear();
}
