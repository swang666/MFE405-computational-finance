#ifndef BINOM_H
#define BINOM_H
#include <stdint.h>

double** question1();
double question2();
double** question3();
double** question4();
double binom_method2(double S0, double K, double u, double d, double p, int n, double T, double r, int type, int put);
double* JR_Btree(double r, double sigma, double dt);
#endif //BINOM_H