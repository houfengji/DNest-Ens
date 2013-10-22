#include <cstdlib>
#include <vector>
#include "max.h"

double max (const std::vector<double> & X) {
	double m = X[0];
	for (size_t i = 1; i < X.size(); ++i) {
		if (X[i] > m) {
			m = X[i];
		}
	}
	return m;
}
