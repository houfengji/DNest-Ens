/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  This header file declares the class data and its members.
 *  This class will have all the data and be used by other classes.
 *
 */

#ifndef DATA_H
#define DATA_H

#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

class Data {
public:
	const string data_name;  // name of the data file
	
	// The following vector are used to store data.
	vector< vector<double> > data;
  
  // The following variables are related to the data
	int64_t data_size;   // size of the data, a.k.a. number of rows
	int64_t num_row;     // just another way to say data_size
	int64_t num_col;     // number of columns
	int64_t k_cutoff;    // only up to k_cutoff rows of data will be used

	// The following members are functions
	int64_t num_of_observations();
	int64_t collect_data();
	
	Data(string);   // string for data file name
	Data(const Data& original);
	~Data() {};
};

#endif
