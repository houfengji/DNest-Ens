/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  In this file, the RV data fitting model is stored, including various kinds
 *  of priors, likelihood, and other related. (J for jitter, D for distiction)
 */
 
#ifndef EXOPLANETJD_H
#define EXOPLANETJD_H


#include <cstdlib>
#include <string>
#include <vector>
#include "data.h"
#include "exoplanet_hyperpara.h"
#include "level.h"
#include "model.h"

class ExoplanetJD: public Model {
public:
	ExoplanetJD(Data &, size_t num_companion, std::string weight_type);
	size_t num_comp;      // number of companions
	size_t num_d;         // number of distinctions, data sources, types equipments...
	double LnPermutation; // the factorial of num_comp
	size_t count_distinction();
	ExoplanetHyperpara hyper;  // hyper-parameters used in priors
	virtual double reparametrize(const std::vector<double> &, std::vector<double> &) const;
	virtual double LnDensity(const std::vector<double> &);
	virtual double LnDensity(const std::vector<double> &, const size_t);
	double LnPrior01(const std::vector<double> &); // all uniform priors
	double LnPrior02(const std::vector<double> &); // more sensible priors
	double LnPrior(const std::vector<double> &, const size_t);
	virtual double LnLikelihood(const std::vector<double> &) ;
	virtual bool orbit_cross (const std::vector<double> &);
	virtual bool physical (const std::vector<double> &);
	virtual bool exchange (const std::vector<double> &);
	void init(size_t ens_size, std::vector< std::vector<double> > & ensemble, double ini);
	virtual ~ExoplanetJD();
private:
	Data data;
};

#endif
