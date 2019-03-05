#ifndef RNG_H
#define RNG_H
#include <stdint.h>

// TODO: add headers that you want to pre-compile here
double * getLGM(int n, int64_t seed);
double * box_muller(double *nums, int n);
double calc_mean(double *nums, int n);
double* vector_minus(double *nums, int n, double num);
double* vector_minus_v(double *nums1, int n, double *nums2, bool add);
double* vector_mult(double *nums, int n, double num);
double vasicek_model_zero_bond(double r0, double sigma, double K, double r_bar, double V, double T, double t, int64_t seed);
double trapzoid_method(int a, int b, double* integrand, double dt);
double cir_model_zero_bond(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t, int64_t seed);
void cir_explicit(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t);
#endif //RNG_H