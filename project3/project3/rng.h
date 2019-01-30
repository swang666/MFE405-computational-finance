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
double* question1(int64_t seed);
double * generate_bivariate(double * z1, double * z2, int n, double rho);
double* question2(int64_t seed);
double * wiener_process(double* nums, int n, double t);
void question3(int64_t seed, double S0, double T, double X, double r, double sigma, double* prices, double* prices_2, double* parta, double** partc, double** partc_2);
double cov(double* X, double* Y, int n);
double* euro_call(double r, double sigma, double S0, double T, double X, double *nums, int n, bool antithetic);
double* question4(int64_t seed);
double normalCDF(double x);
double black_schole(double r, double sigma, double S0, double T, double X);
double* geometric_brownian_motion(double r, double sigma, double S0, double T, double* nums, int n);
void question5(int64_t seed, double* parta, double** partb, double* partc, double** partd);
double* brownian_motion_path(double r, double sigma, double S0, double T, double* nums, int n);
double euler_approx(double t, int n, double x0);
double myfunc(double x);
double* question6(int64_t seed);
double* monte_carlo(double* nums, int n, double(*f)(double));
double monte_carlo_euro_call(double S0, double T, double X, double r, double sigma, int64_t seed);
double sim_normalCDF(double x);
double pdf(double x);
double* greeks(double S0, double T, double r, double X, double sigma);
double* greeks_2(double S0, double T, double r, double X, double sigma, int64_t seed);
double heston_model(double rho, double r, double S0, double V0, double sigma, double alpha, double beta, double T, double K, int type, int64_t seed);
#endif //RNG_H