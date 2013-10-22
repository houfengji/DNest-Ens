#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include "mean.h"

using namespace std;

double mean (const vector<double> & X) {
	double m = 0;
	for (size_t i = 0; i < X.size(); ++i) {
		m += X[i];
	}
	m /= static_cast<double>(X.size());
	return m;
}

double logMeanLog (const vector<double> & logX) {
	double maxLogX = *max_element(logX.begin(), logX.end());
	double lm = 0;
	for (size_t i = 0; i < logX.size(); ++i) {
		lm += exp(logX[i] - maxLogX);
	}
	lm /= static_cast<double>(logX.size());
	lm = log(lm) + maxLogX;
	return lm;
}

