#ifndef MBS_H
#define MBS_H
#include <stdint.h>

double * getLGM(int n, int64_t seed);
double * box_muller(double *nums, int n);
double calc_mean(double *nums, int n);
double trapzoid_method(int a, int b, double* integrand, double dt);
double MBS(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double x, int64_t seed);
double OAS(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double mkt_p, int64_t seed);
double MBS_explicit(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double x);
double OAS_explicit(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double mkt_p);
#endif //MBS_H

