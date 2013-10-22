/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *
 *  In this cpp file, all the member functions of class Data are defined.
 *
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include "data.h"

using namespace std;


// This function checks how many lines there are in the data file and returns
// the number of lines which is also the number of observations. The return
// value is recorded in class member data_size.
int64_t Data::num_of_observations() {
	string data_file_name("./data/" + data_name + ".dat"); 
	ifstream input(data_file_name.c_str());

	if(input.fail()) {
		cerr << "error: unable to open input data file: " << data_file_name << endl;
		return -1;
	}

	int64_t num_of_lines = 0;
	string temporary; // temporary string used to store the line read by getline()
	while ( !input.eof() ) {
		getline( input, temporary ); // get and count a line until the line is EOF
		if( temporary[0] != 0 ) ++ num_of_lines;
	}

	input.close();
	return num_of_lines;
}


// This function reads the data into 2-d vector 'data'.
int64_t Data::collect_data() {
	string data_file_name("./data/" + data_name + ".dat"); 
	ifstream input(data_file_name.c_str());
  
	if(input.fail()) {
		cerr << "error: unable to open input data file: " << data_file_name << endl;
		return -1;
	}
	
	vector<double> input_array; // used to store the data in an array
 	double temp;
 	while(!input.eof()) {
 		input >> temp;
		input_array.push_back(temp);
	}
	
	num_col = input_array.size() / data_size; // determining number of columns
	data.resize(data_size, vector<double>(num_col, 0));
	for (int64_t i = 0; i < data_size; ++i) {
		for (int64_t j = 0; j < num_col; ++j) {
			data[i][j] = input_array[i*num_col + j];
		}
	}
  
	input.close();
	return 1;
}

// constructor
Data::Data(string name):          // name of the data
           data_name(name)
{

	data_size = num_of_observations();
	num_row = data_size;
	if (data_size > 0) {
		collect_data();
	}
	else {
		cerr << "No Data Found!" << endl;
	}
	
}

Data::Data(const Data& original):
           data_name(original.data_name),
           num_col(original.num_col),
           num_row(original.num_row),
           k_cutoff(original.k_cutoff),
           data_size(original.data_size)
{
	if (data_size > 0) {
		data.resize(data_size, vector<double>(num_col, 0));
		for (int64_t i = 0; i < data_size; ++i) {
			for (int64_t j = 0; j < num_col; ++j) {
				data[i][j] = original.data[i][j];
			}
		}
	}
}
