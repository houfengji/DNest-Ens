/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  This header file declares a model with uniform prior.
 *  The purpose of this model is to check the performance of the algorithm.
 *
 */
 
#ifndef UNIFORMTEST_H
#define UNIFORMTEST_H

#include <cstdlib>
#include <string>
#include <vector>
#include "data.h"
#include "level.h"
#include "model.h"

class UniformTest: public Model {
public:
	UniformTest(Data &, size_t dimension, std::string weight_type);
	double lower_bound;  // hyper-parameter
	double upper_bound;  // hyper-parameter
	virtual double reparametrize(const std::vector<double> &, std::vector<double> &) const;
	virtual double LnDensity(const std::vector<double> &) ;
	virtual double LnDensity(const std::vector<double> &, const size_t) ;
	double LnPrior(const std::vector<double> &) ;
	double LnPrior(const std::vector<double> &, const size_t) ;
	virtual double LnLikelihood(const std::vector<double> &);
	void init(size_t, std::vector< std::vector<double> > &, double) const;
	double theo_lnprims_2d(double ln_likelihood);
	virtual ~UniformTest();
private:
	Data data;
};

#endif
