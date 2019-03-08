#ifndef RNG_H
#define RNG_H
#include <stdint.h>

// TODO: add headers that you want to pre-compile here
double * getLGM(int n, int64_t seed);
double * box_muller(double *nums, int n);
double calc_mean(double *nums, int n);
double vasicek_model_zero_bond(double r0, double sigma, double K, double r_bar, double V, double T, double t, int64_t seed);
double vasicek_model_semicoupon_bond(double r0, double sigma, double C, double K, double r_bar, double V, double T, double t, int64_t seed);
double vasicek_model_zero_bond_call(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t, int64_t seed);
double vasicek_model_coupon_bond_call(double r0, double sigma, double C, double K, double strike, double r_bar, double V, double T, double S, double t, int64_t seed);
double trapzoid_method(int a, int b, double* integrand, double dt);
double cir_model_zero_bond(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t, int64_t seed);

double cir_explicit(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t);
double G2pp(double r0, double x0, double y0, double phi0, double rho, double a, double b, double sigma, double eta, double phi_t, double strike, double V, double T, double S, double t);
double G2pp_mc(double r0, double x0, double y0, double phi0, double rho, double a, double b, double sigma, double eta, double phi_t, double strike, double V, double T, double S, double t, int64_t seed);
double * generate_bivariate(double * z1, double * z2, int n, double rho);
double normalCDF(double x);
#endif //RNG_H