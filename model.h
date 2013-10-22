/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  Jan 27, 2013
 *  This header file declares the abstract base class Model.
 *  All models will be derivatives of this class.
 *  This class has some basic information about the model, i.e.
 *  dimension, prior, likelihood etc.
 *
 */
 
#ifndef MODEL_H
#define MODEL_H

#include <cstdlib>
#include <string>
#include <vector>
#include "level.h"

class Model {
public:
	Model(size_t dimension,          // the dimension of the model's paramter space
	      std::string weight_type);  // weight type for the level
	Model(std::string weight_type);  // weight type for the level
	size_t dim;                      // dimension of the model's parameter space
	std::vector<double> best_fit;    // This vector is used to keep track of the best fit parameters
	Level level;                     // constructor of level requires 'weight type'
	std::string model_name;          // name of the model, for output file names
	std::string time_label;          // the time when an object is made, used to label output files
	std::vector<double> chain_LnLikelihood;   // kuh mjieng si ngje
	size_t total_num_visit;
	std::vector<size_t> Level_Visits_T;       // The visits to each level including burning phase
	std::vector<size_t> Level_Visits;         // The visits to each level
	std::vector<size_t> Above_Visits;         // The visits of the level above while visiting each level
	std::vector<double> Quant_Ens_Mean;       // above/visit mean of the ensemble, used to evaluate tau
	                                          // Threshold k+1, Quantizations[k]=1; otherwise, Quantizations[k]=0.
	double evidence;
	double evidence_r;               // refined evidence
	double evidence_err1;            // error bar by integrand
	double evidence_err2;            // error bar by prior mass
	double evidence_err3;            // error bar by prior mass cross terms
	virtual double reparametrize(const std::vector<double> &, std::vector<double> &) const = 0;
	virtual double LnDensity(const std::vector<double> &) = 0;
	virtual double LnDensity(const std::vector<double> &, const size_t) = 0;
	virtual double LnLikelihood(const std::vector<double> &) = 0;
	void clear_all_visits(void);
	void clear_for_sampling(void);
	virtual ~Model() {}
private:
};

#endif
