#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include "model.h"

double sampling (Model & model,      // data and model
                 std::vector< std::vector<double> > & ensemble,
                 const size_t num_step,
                 const double a);

#endif
