/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  Jan 27, 2013
 *  This cpp file defines all the members of Level class
 *
 */
 
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "exception.h"
#include "level.h"

using namespace std;

Level::Level(const std::string & type): weight_type(type)
{
	num_level = 1;      // start with one level M_0 = M_0^* = 1, L_0^* = 0
	LWeight_norm = 0;   // ensure sum of weights is 1
	exponential_weight_coefficient = 0.5; // go down 2 levels with Prob = exp(-1)
	LnThres.push_back(-numeric_limits<double>::infinity()); // L_0^* = 0
	LnPrims.push_back(0.);                                  // M_0 = M_0^* = 1
	bool updated = update_lweight();
	if (updated == 0) {
		throw ( Exception("Level, Constructor : Updating level weights failed!") );
	}
}

void Level::add_level(double L, double X) {
	++ num_level;
	LnThres.push_back(L);
	LnPrims.push_back(X);
	bool updated = update_lweight();
	if (updated == 0) {
		throw ( Exception("Level, add_level : Updating level weights failed!") );
	}
}

void Level::change_weight_type (const std::string & type) {
	bool updated;
	weight_type = type;
	updated = update_lweight();
	if (updated == 0) {
		throw ( Exception("Level, change_weight_type : Updating level weights failed! Type " + type) );
	}
}

void Level::clear_all() {
	num_level = 1;
	LnThres.resize(1);
	LnPrims.resize(1);
	LnThres[0] = -numeric_limits<double>::infinity();
	LnPrims[0] = 0.;
	bool updated = update_lweight();
	if (updated == 0) {
		throw ( Exception("Level, clear_all : Updating level weights failed!") );
	}
}

void Level::pop_one_level() {
	num_level = num_level - 1;
	LnThres.pop_back();
	LnPrims.pop_back();
	bool updated = update_lweight();
	if (updated == 0) {
		throw ( Exception("Level, pop_one_level : Updating level weights failed!") );
	}
}

// Returned value is whether weight updating succeeded.
bool Level::update_lweight() {
	LWeight.resize(LnThres.size());
	
	if (weight_type == "exponential") {
		double J = (double)LnThres.size() - 1;
		double  Lambda = exponential_weight_coefficient;
		for (size_t j = 0; j < LnThres.size(); ++j) {
			LWeight[j] = ((double)j-J)*Lambda;
		}
		LWeight_norm = 0;
		for (size_t i = 0; i < LnThres.size(); ++i) {
			LWeight_norm += exp(LWeight[i]);
		}
		LWeight_norm = log(LWeight_norm);
		return 1;
	}
	if (weight_type == "delta") {
		for (int64_t j = 0; j < LnThres.size()-1; ++j) {
			LWeight[j] = -1e308;
		}
		LWeight[LnThres.size()-1] = 0;
		LWeight_norm = 0;
		return 1;
	}
	if (weight_type == "uniform") {
		for (int64_t j = 0; j < LnThres.size(); ++j) {
			LWeight[j] = 0;
		}
		LWeight_norm = 0;
		for (size_t i = 0; i < LnThres.size(); ++i) {
			LWeight_norm += exp(LWeight[i]);
		}
		LWeight_norm = log(LWeight_norm);
		return 1;
	}
	std::cerr << "Weight type doesn't exist!" << std::endl;
	return 0;
}

