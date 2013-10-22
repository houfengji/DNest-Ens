/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  Jan 27, 2013
 *  This cpp file defines a simple constructor for Level class.
 *
 */
 
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "level.h"

Level::Level(const std::string & type): weight_type(type)
{
	num_level = 1;
	LWeight_norm = 0;
	LnThres.push_back(-1.0e300);     // no threshold at all
	LnPrims.push_back(0.);           // The prior mass the lowest level covers is 1.
	bool updated = update_lweight();
	exponential_weight_coefficient = 0.5;
	if (updated == 0) {
		std::cerr << "Level, Constructor : Updating level weights failed!" << std::endl;
	}
}

void Level::add_level(double L, double X) {
	++ num_level;
	LnThres.push_back(L);
	LnPrims.push_back(X);
	bool updated = update_lweight();
	if (updated == 0) {
		std::cerr << "Level, add_level : Updating level weights failed!" << std::endl;
	}
}

void Level::change_weight_type (const std::string & type) {
	bool updated;
	weight_type = type;
	updated = update_lweight();
	
	if (updated == 0) {
		std::cerr << "Level, change_weight_type : Updating level weights failed! Type " << type << std::endl;
	}
}

void Level::clear_all() {
	num_level = 1;
	LnThres.resize(1);
	LnPrims.resize(1);
	bool updated = update_lweight();
	if (updated == 0) {
		std::cerr << "Level, Constructor : Updating level weights failed!" << std::endl;
	}
}

void Level::pop_one_level() {
	if (LnLikelihood_var.size() == num_level - 1) {
		LnLikelihood_var.pop_back();
	}
	num_level = num_level - 1;
	LnThres.pop_back();
	LnPrims.pop_back();
	bool updated = update_lweight();
	if (updated == 0) {
		std::cerr << "Level, Constructor : Updating level weights failed!" << std::endl;
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

