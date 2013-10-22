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
