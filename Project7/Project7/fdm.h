#ifndef FDM_H
#define FDM_H
#include <stdint.h>
#include <Eigen>

using namespace Eigen;

MatrixXd EFDM(double S0, int k, int euro, int call);
MatrixXd IFDM(double S0, int k, int euro, int call);
MatrixXd CNFDM(double S0, int k);
MatrixXd GFD(double S0, double dS, double alpha, int call);
#endif //FDM_H