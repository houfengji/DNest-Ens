/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  The class declared in this header is used to store best-fit parameters, or
 *  posterior-median parameters, to initialize the ensemble. They can be as
 *  specific as one wants, but don't have to be.
 */

#ifndef EXOPLANET_INIT_H
#define EXOPLANET_INIT_H

#include <string>
#include <vector>

struct ExoplanetInit {
	ExoplanetInit(std::string data_name, size_t num_comp, size_t num_d, size_t id);
	ExoplanetInit(std::string data_name, size_t num_comp, size_t num_d);
	ExoplanetInit(std::string data_name, size_t num_comp);
	
	std::vector<double> init;
};


#endif
