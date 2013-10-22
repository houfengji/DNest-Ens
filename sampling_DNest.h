#ifndef SAMPLING_DNEST_H
#define SAMPLING_DNEST_H

#include "model.h"
#include "preprocessors.h"

double sampling_DNest (Model & model,
                       const size_t num_step,
                       const size_t step_size,
                       std::vector< std::vector<double> > & ensemble,
                       std::vector<size_t> & levels,
                       double a);

bool new_level (Model & model,
                std::vector< std::vector<double> > & ensemble,
                std::vector<size_t> & levels,
                const double a,
                const size_t burn_in_steps,   // number of steps before actually building the chain
                const size_t chain_length,    // length of the chain used to determine next level
                const size_t chain_length_,    // length of the chain used to determine next level
                const size_t step_size);      // number of steps between two entries in the chain

double LLmean(Model & model, matrix & ensemble);
void ensemble_mean(const matrix & ensemble, std::vector<double> & mean);

double accept_ensemble (Model & model, 
                        const size_t & level,
                        const matrix & ensemble, 
                        const matrix & proposed_ensemble);
                        
size_t update_ensemble (Model & model,      // data and model
                        std::vector< std::vector<double> > & ensemble,
                        const std::vector<size_t> & levels,
                        const double a);

size_t update_index (Model & model,
                     const std::vector< std::vector<double> > & ensemble,
                     std::vector<size_t> & levels);
                  
#endif
