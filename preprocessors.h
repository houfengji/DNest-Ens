#ifndef PREPROCESSORS_H
#define PREPROCESSORS_H

#define TESTING 999
// TESTING = 1 for outputting all the samples when building levels, so able to see where the thresholds and other properties of the levels.
// TESTING = 2 for outputting all the samples to a file after all levels have been built.
// TESTING = 3 for outputting all the likelihoods and corresponding X-Values after assigning x values.
// TESTING = 4 for outputting all the likelihood samples after levels have been built but keep the likelihoods that are between two levels in separate files.

#include <vector>

typedef std::vector< std::vector<double> > matrix;

#endif
