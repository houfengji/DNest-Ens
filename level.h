/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  Oct 25, 2013
 *  This header includes the declaration of Level class. This class stores all
 *  the levels' information while constructing the levels and also the refined
 *  levels. This class has the weight of each level as well.
 *
 */


#ifndef LEVEL_H
#define LEVEL_H

#include <cstdlib>
#include <string>
#include <vector>

class Level {
public:
	Level(const std::string &);    // The string argument is the type of weight.
	std::vector<double>::size_type num_level;
	std::string weight_type;
	std::vector<double> LnThres;             // likelihood thresholds L_j^*
	std::vector<double> LnPrims;             // prior masses M_j
	std::vector<double> LWeight;             // weight of levels 
	std::vector<double> LnPrims_r;           // refined prior masses M_j^*
	std::vector<double> LnPrims_var;         // variance of M_j^*
	std::vector<double> LnInterval_var;      // variance of samples between levels
	double LWeight_norm;                     // normalization of the weights
	double exponential_weight_coefficient;   // only used in exponential weight
	void add_level(double L, double X);
	void clear_all();
	void pop_one_level();
	bool update_lweight();
	void change_weight_type(const std::string &);
	~Level() {}
private:
};

#endif
