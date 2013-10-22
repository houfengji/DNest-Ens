/*
 *   Fengji Hou
 *   fh417@nyu.edu
 *   New York University
 *   This cpp file contains a simple function which takes float/double number
 *   and returns the sign of it (as int).
 *
 */

#include "sign.h"
 
int sign (const double value) {
	
	if (value < 0) {
		return -1;
	}
	else {
		return 1;
	}
}
