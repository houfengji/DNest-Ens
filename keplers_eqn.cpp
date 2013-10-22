/*
 *   Fengji Hou
 *   fh417@nyu.edu
 *   New York University
 *   This cpp file contains the Kepler's equation, the orbital equations
 *   for the orbital parameters.
 *
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "keplers_eqn.h"

using namespace std;

void eccentric_anomaly_test (double M, double e, stringstream & label_ss);

int sign (const double value) {
	
	if (value < 0) {
		return -1;
	}
	else {
		return 1;
	}
}

double mean_anomaly (const double E,       //eccentric anomaly
                     const double e) {     //eccentricity
	
	return E - e * sin(E);
}


double eccentric_anomaly (const double rM,     //mean anomaly
                          const double e) {    //eccentricity
	
	double M = rM - floor(rM/(2.0*M_PI))*2*M_PI;
	double k = floor(M/M_PI); // If 0<M<=pi, which = 0. If pi<M<=2pi, which = 1;
	double E = M + e * sin(M);
	double Emin = k*M_PI;
	double Emax = (k+1)*M_PI;
	double Ep;
	size_t iteration = 0;
	double delta_M = 1.0;
	
	while (fabs(delta_M/M) > 0.00000001) {
		delta_M = M - (E - e * sin(E));
		Ep = E + delta_M / (1.0 - e * cos(E));
		++iteration;
		if (Ep < Emin) {
			E = 0.5 * (E+Emin);
		}
		else {
			if (Ep > Emax) {
				E = 0.5 * (E+Emax);
			}
			else {
				E = Ep;
			}
		}
	}
	return E;
}

double true_anomaly (const double E,        //eccentric anomaly
                     const double e) {      //eccentricity
	
	double f;
	f = acos((cos(E) - e) / (1.0 - e * cos(E)));
	f *= sign(sin(f)) * sign(sin(E)); // sin(f) will always have the same sign as sin(E) does.
	return f;
}

// radial_velocity returns the radial velocity of the star on our line of view,
// if we know all the companion orbit data.
double rad_v (const double A,          //amplitude
              const double f,          //true anomaly
              const double e,          //eccentricity
              const double cpi) {      //2nd phase
	
	return A * (sin(f + cpi) + e * sin(cpi));
}

// radial_velocity_predicted returns the radial velocity based on the machine learning parameters
// The arguments of this function are the 5 orbital parameters for MCMC
double rad_v_pred (const double t,           //time
                   const double amplitude,
                   const double omega,       //angular speed
                   const double phi,         //1st phase
                   const double e,           //eccentricity
                   const double cpi) {       //2nd phase
	
	double M = omega * t + phi;
	double E = eccentric_anomaly(M, e);
	double f = true_anomaly(E, e);
	double rv = rad_v(amplitude, f, e, cpi);
	//if (isnan(rv)) {
		//rad_v_test(t, amplitude, omega, phi, e, cpi);
	//}
	return rv;
}

void rad_v_test (const double t,           //time
                 const double amplitude,
                 const double omega,       //angular speed
                 const double phi,         //1st phase
                 const double e,           //eccentricity
                 const double cpi) {       //2nd phase
	
	fstream out;
	stringstream label_ss (stringstream::in | stringstream::out);
	label_ss << time(NULL);
	label_ss << "_";
	label_ss << rand();
	out.open(("Kepler_Eqn_Exception_"+label_ss.str()+".txt").c_str(), ios::out);
	out << "time         = " << setprecision(16) << t << endl;
	out << "amplitude    = " << setprecision(16) << amplitude << endl;
	out << "omega        = " << setprecision(16) << omega << endl;
	out << "phi          = " << setprecision(16) << phi << endl;
	out << "eccentricity = " << setprecision(16) << e << endl;
	out << "varpi        = " << setprecision(16) << cpi << endl;
	double M = omega * t + phi;
	out << "mean anomaly = " << setprecision(16) << M << endl;
	eccentric_anomaly_test(M, e, label_ss);
}

void eccentric_anomaly_test (double M, double e, stringstream & label_ss) {
	double E = M + e * sin(M);
	int64_t iteration = 0;
	double delta_M = 1.0;
	fstream out;
	out.open(("Ecc_Anomaly_Exception_"+label_ss.str()+".txt").c_str(), ios::out);
	out << setprecision(16) << "Mean Anomaly = " << M << "   Ecc = " << e << endl;
	out << "  E    cos(E)   e*cos(E)   delta_M" << endl;
	//Because M = E - e * sin(E) cannot be solved analytically, we use iteration.
	while (fabs(delta_M / M) > 0.000001) { 
		delta_M = M - mean_anomaly(E, e);
		out << setprecision(16) << E << "   " << cos(E) << "   " << e * cos(E) << "   " << delta_M << endl;
		E = E + delta_M / (1.0 - e * cos(E));
		iteration += 1;
		if (iteration > 100000 && delta_M < 0.1) {
			break;
			//Sometimes, iteration can get stuck mostly due to high eccentricity.
			//This section helps the code jump out of the loop
		}
	}
}
