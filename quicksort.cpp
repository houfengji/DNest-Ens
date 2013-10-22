/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  This cpp file contains quick-sort sorting algorithm.
 *  Sorted array will be in descending order.
 *
 */

#include <iostream>
#include <vector> 

void quicksort(std::vector<double> & vec, 
               std::vector<double>::iterator left,
               std::vector<double>::iterator right) {

	std::vector<double>::iterator i = left;
	std::vector<double>::iterator j = right;
	std::vector<double>::difference_type size = j - i;
	double pivot = *(i + size/2);
	
	double tmp;
	
	while (i <= j) {
		while (*i > pivot) {
			++i;
		}
		while (*j < pivot) {
			--j;
		}
		if (i <= j) {
			tmp = *i;
			*i  = *j;
			*j  = tmp;
			++i;
			--j;
		}
	}

	if (left < j) {
		quicksort (vec, left, j);
	}
	if (i < right) {
		quicksort (vec, i, right);
	}
} 
