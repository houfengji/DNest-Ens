#ifndef ACOR_H
#define ACOR_H

#include <cstdlib>
#include <vector>

int acor( double & mean, double & sigma, double & tau, std::vector<double> X, size_t);
int logacorlog( double & logmean, double & logsigma, double & tau, std::vector<double> & logX);

#endif
