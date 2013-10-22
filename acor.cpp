#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include "acor.h"
#include "exception.h"
#include "int2str.h"

using namespace std;

#define TAUMAX  100              //  Compute tau directly only if tau < TAUMAX. Otherwise compute tau using the pairwise sum series.
#define WINMULT 10               //  Compute autocovariances up to lag s = WINMULT*TAU
#define MAXLAG  TAUMAX*WINMULT   //  The autocovariance array is double C[MAXLAG+1] so that C[s] makes sense for s = MAXLAG.
#define MINFAC   5               //  Stop and print an error message if the array is shorter than MINFAC * MAXLAG.

/*  Jonathan Goodman, March 2009, goodman@cims.nyu.edu  */
/*  Adapted by Fengji Hou, Feb 2013, fh417@nyu.edu, again on Jul 2013. */

int acor( double & mean, double & sigma, double & tau, vector<double> X, size_t L){
	
	static int recursive_depth = 0;
	++recursive_depth;
	mean = 0.;                                   // Compute the mean of X ... 
	for ( size_t i = 0; i < L; ++i) {
		mean += X[i];
	}
	mean /= static_cast<double>(L);
	
	for ( size_t i = 0; i < L; ++i ) {
		X[i] -= mean;                              //  ... and subtract it away.
	}
   
	if ( L < MINFAC*MAXLAG ) {
		recursive_depth = 0;
		throw( Exception("acor: Too Many Lags Are Required to Evaluate Autocovariance!") );
	}
	
	double C[MAXLAG+1]; 
	for ( int s = 0; s <= MAXLAG; ++s ) {
		C[s] = 0.;                                // Here, s=0 is the variance, s = MAXLAG is the last one computed.
	}
	
	size_t iMax = L - MAXLAG;                   // Compute the autocovariance function . . . 
	for ( size_t i = 0; i < iMax; ++i ) {

		for ( int s = 0; s <= MAXLAG; ++s ) {
			C[s] += X[i]*X[i+s];                    // ...  first the inner products ...
		}
	}
	
	for ( int s = 0; s <= MAXLAG; ++s ) {
		C[s] = C[s]/static_cast<double>(iMax);    // ...  then the normalization.
		//cout << s << "    " << C[s] << endl;
	}

	double D = C[0];                            // The "diffusion coefficient" is the sum of the autocovariances
	for ( int s = 1; s <= MAXLAG; ++s ) {
		D += 2*C[s];                              // The rest of the C[s] are double counted since C[-s] = C[s].
	}
	
	if (D < C[0]) {
		
		D = C[0];
	}
	sigma = sqrt( D / static_cast<double>(L) ); // The standard error bar formula, if D were the complete sum.
	tau   = D / C[0];                           // A provisional estimate, since D is only part of the complete sum.
  

	if ( tau*WINMULT < MAXLAG ) {
		cout << "Recursive Depth = " << recursive_depth << endl;
		recursive_depth = 0;
		return 0;             // Stop if the D sum includes the given multiple of tau.
                          // This is the self consistent window approach.
	}
	else {                  // If the provisional tau is so large that we don't think tau
                          // is accurate, apply the acor procedure to the pairwase sums
											    // of X.
		size_t Lh = L/2;                                // The pairwise sequence is half the length (if L is even)
		double newMean;                                 // The mean of the new sequence, to throw away.
		int j1 = 0;
		int j2 = 1;
		for ( size_t i = 0; i < Lh; ++i ) {
			X[i] = X[j1] + X[j2];
			j1  += 2;
			j2  += 2; 
		}
	  acor( newMean, sigma, tau, X, Lh);
	  D      = .25*(sigma) * (sigma) * L;             // Reconstruct the fine time series numbers from the coarse series numbers.
	  tau   = D/C[0];                                 // As before, but with a corrected D.
	  sigma = sqrt( D/ static_cast<double>(L));       // As before, again.
	}

	return 0;
}

int logacorlog( double & logmean, double & logsigma, double & tau, vector<double> & logX) {
	size_t L = logX.size();
	double max_lX = *(max_element(logX.begin(), logX.end()));
	vector<double> X(L,0);
	for (size_t i = 0; i < L; ++i) {
		X[i] = exp(logX[i] - max_lX);
	}
	double mean, sigma;
	acor(mean, sigma, tau, X, L);
	logsigma = log(sigma) + max_lX;
	logmean  = log(mean)  + max_lX;
	return 1;
}
