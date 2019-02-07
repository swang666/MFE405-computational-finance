#ifndef BINOM_H
#define BINOM_H
#include <stdint.h>

double** question1();
double question2();
double** question3();
double** question4();
double** question5();
double question6();

double binom_method2(double S0, double K, double u, double d, double p, int n, double T, double r, int type, int put);
double* JR_Btree(double r, double sigma, double dt);
double trinomial_method(double S0, double K, double u, double d, double pu, double pd, double pm, int n, double T, double r);
double trinomial_method2(double X0, double K, double delta_Xu, double delta_Xd, double pu, double pd, double pm, int n, double T, double r);
double* halton_seq(int base, int num);
double * box_muller(double *nums1, double* nums2, int n);
double halton_euro_call(double S0, double K, double T, double r, double sigma, int n, int b1, int b2);
double calc_mean(double *nums, int n);
#endif //BINOM_H