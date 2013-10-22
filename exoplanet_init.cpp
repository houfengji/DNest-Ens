#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "exoplanet_init.h"

using namespace std;

ExoplanetInit::ExoplanetInit(string data_name, size_t num_comp) {
	
	bool found = 0;
	
	if (data_name.compare("x001") == 0 && num_comp == 0) {
		found = 1;
		init.resize(num_comp * 5 + 2);
		init[0] = 0.3956;
		init[1] = 3532.13;
	}
	
	if (data_name.compare("x001") == 0 && num_comp == 1) {
		found = 1;
		init.resize(num_comp * 5 + 2);
		init[0] = 100;
		init[1] = 0.12566370614;
		init[2] = 2;
		init[3] = 0.1;
		init[4] = 4;
		init[5] = 5;
		init[6] = 25;
	}
	
	if (data_name.compare("x001") == 0 && num_comp == 2) {
		found = 1;
		init.resize(num_comp * 5 + 2);
		init[0] = 8.72;
		init[1] = 0.2950519232;
		init[2] = 2.87;
		init[3] = 0.051;
		init[4] = 4.81;
		init[5] = 100;
		init[6] = 0.12566370614;
		init[7] = 1.82;
		init[8] = 0.092;
		init[9] = 4.15;
		init[10] = 5.8;
		init[11] = 0.00052;
	}
	
	if (found == 0) {
		cerr << "Exoplanet Initialization Failed!" << endl;
		cerr << "Data Name or Model couldn't be found!" << endl;
	}
}

ExoplanetInit::ExoplanetInit(string data_name, size_t num_comp, size_t num_d) {
	
	bool found = 0;
	
	if (found == 0) {
		cerr << "Exoplanet Initialization Failed!" << endl;
		cerr << "Data Name or Model couldn't be found!" << endl;
	}
}
