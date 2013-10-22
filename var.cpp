#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include "mean.h"
#include "var.h"

using namespace std;

double var(const vector<double> & X) {
	double m = mean(X);
	double v = 0;
	for (size_t i = 0; i < X.size(); ++i) {
		v += X[i] * X[i];
	}
	v = v / static_cast<double>(X.size()) - m*m;
	return v;
}

double logVarLog(const vector<double> & logX) {
	double maxLogX = *max_element(logX.begin(), logX.end());
	double lm = logMeanLog(logX) - maxLogX;
	double lv = 0;
	for (size_t i = 0; i < logX.size(); ++i) {
		lv += exp( 2.0 * (logX[i] - maxLogX) );
	}
	lv = lv / static_cast<double>(logX.size()) - exp(2.0 * lm);
	lv = log(lv) + 2.0 * maxLogX;
	return lv;
}

