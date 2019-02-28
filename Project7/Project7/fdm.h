#ifndef FDM_H
#define FDM_H
#include <stdint.h>


double EFDM(double S0, int k, int euro, int call);
double IFDM(double S0, int k, int euro, int call);
double CNFDM(double S0, int k);
double EFDM_S(double S0, double dS, double alpha, int call);
#endif //FDM_H