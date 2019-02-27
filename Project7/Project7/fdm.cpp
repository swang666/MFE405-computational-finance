#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "fdm.h"
#include "matrix.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;

double EFDM(double S0, int k) {
	double sigma = 0.2;
	double T = 0.5;
	double K = 10;
	double dt = 0.002;
	double N = T / dt;
	double dx = sigma * sqrt(k*dt);
	double r = 0.04;
	double Smax = 16;
	double Smin = 4;
	double pu = dt * (sigma*sigma / (2 * dx*dx) + (r - sigma * sigma / 2) / (2 * dx));
	double pm = 1 - dt * sigma*sigma / (dx*dx) - r * dt;
	double pd = dt * (sigma*sigma / (2 * dx*dx) - (r - sigma * sigma / 2) / (2 * dx));
	int num_path = (log(S0) - log(Smin)) / dx;
	double** XN = zero_matrix(2 * num_path + 1, 1);
	double** FN = zero_matrix(2 * num_path + 1, 1);
	for (int i = num_path; i < 2 *num_path + 1; ++i) {
		XN[i][0] = log(S0) - (i - num_path) * dx;
		FN[i][0] = max(K - exp(XN[i][0]), 0.0);
	}
	for (int i = num_path-1; i >= 0; --i) {
		XN[i][0] = log(S0) + (num_path - i) * dx;
		FN[i][0] = max(K - exp(XN[i][0]), 0.0);
	}
	int total_num_path = 2 * num_path + 1;
	double** A = zero_matrix(total_num_path, total_num_path);
	A[0][0] = pu;
	A[0][1] = pm;
	A[0][2] = pd;
	A[total_num_path-1][total_num_path-3] = pu;
	A[total_num_path-1][total_num_path-2] = pm;
	A[total_num_path-1][total_num_path-1] = pd;
	for (int i = 1; i < total_num_path-1; ++i) {
		A[i][i - 1] = pu;
		A[i][i] = pm;
		A[i][i + 1] = pd;
	}

	double** Fi;
	double** temp;
	Fi = matrix_add(zero_matrix(total_num_path, 1), FN, total_num_path, 1);
	double** B = zero_matrix(total_num_path, 1);
	B[total_num_path - 1][0] = -exp(XN[total_num_path - 1][0]) + exp(XN[total_num_path - 2][0]);
	for (int i = 0; i < N; ++i) {
		temp = matrix_add(matrix_mult(A, Fi, total_num_path, total_num_path, total_num_path, 1), B, total_num_path, 1);
		Fi = matrix_add(zero_matrix(total_num_path, 1), temp, total_num_path, 1);
	}
	return Fi[num_path][0];
}