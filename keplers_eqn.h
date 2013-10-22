#ifndef KEPLERS_EQN_H
#define KEPLERS_EQN_H

double mean_anomaly (const double E, const double e);

double eccentric_anomaly (const double M, const double e);

double true_anomaly (const double E, const double e);

double rad_v (const double A, const double f, const double e, const double cpi);

double rad_v_pred (const double t, const double amplitude, const double omega, const double phi, const double e, const double cpi);

void rad_v_test (const double t, const double amplitude, const double omega, const double phi, const double e, const double cpi);

#endif
