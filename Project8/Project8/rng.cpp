#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "rng.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;

const double pi = atan(1) * 4;

double * getLGM(int n, int64_t seed) {
	int64_t m = pow(2, 31) - 1;
	int a = pow(7, 5);
	int b = 0;
	int64_t X0 = seed; //some arbitrary number
	double* r = new double[n];
	int64_t* x = new int64_t[n];
	x[0] = X0;
	r[0] = (double)(X0 + 0.5) / m;
	for (int i = 1; i < n; ++i) {
		x[i] = (a * x[i - 1]) % m;
		r[i] = (double)(x[i] + 0.5) / m;
	}
	return r;
}

double * box_muller(double *nums, int n) {
	double z1, z2;
	double *std_norm = new double[n];
	for (int i = 0; i < n; i = i + 2) {
		z1 = sqrt(-2 * log(nums[i]))* cos(2 * pi * nums[i + 1]);
		z2 = sqrt(-2 * log(nums[i]))* sin(2 * pi * nums[i + 1]);
		std_norm[i] = z1;
		std_norm[i + 1] = z2;
	}
	return std_norm;
}

double calc_mean(double *nums, int n) {
	double sum_of_mean = 0;
	for (int i = 0;i < n;++i) {
		sum_of_mean = sum_of_mean + nums[i];
	}
	double mean = sum_of_mean / n;
	return mean;
}

double* vector_minus(double *nums, int n, double num) {
	double* result = new double[n];
	for (int i = 0;i < n;++i) {
		result[i] = nums[i] - num;
	}
	return result;
}

double* vector_minus_v(double *nums1, int n, double *nums2, bool add) {
	double* result = new double[n];
	for (int i = 0;i < n;++i) {
		if (add == false) {
			result[i] = nums1[i] - nums2[i];
		}
		else {
			result[i] = nums1[i] + nums2[i];
		}
	}
	return result;
}

double* vector_mult(double *nums, int n, double num) {
	double* result = new double[n];
	for (int i = 0; i < n;++i) {
		result[i] = nums[i] * num;
	}
	return result;
}


double vasicek_model_zero_bond(double r0, double sigma, double K, double r_bar, double V, double T, double t, int64_t seed) {

	double* result = new double[10000];
	for (int j = 0; j < 10000; ++j) {
		int n = 126;
		double dt = (T - t) / n;
		double* rs = new double[n + 1];
		rs[0] = r0;
		double* rand_num = getLGM(n, seed+j*100);
		double* z = box_muller(rand_num, n);
		for (int i = 1; i < n + 1; ++i) {
			rs[i] = rs[i - 1] + K * (r_bar - rs[i - 1])*dt + sigma * sqrt(dt)*z[i - 1];
		}
		double P = 1000 * exp(-trapzoid_method(0, n, rs, dt));
		result[j] = P;
		delete[] rand_num;
		delete[] z;
		delete[] rs;
	}
	return calc_mean(result, 10000);
}

double cir_model_zero_bond(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t, int64_t seed) {

	double* result = new double[10000];
	int n = 252;
	double dt = (S - t) / n;
	cout << dt << endl;
	for (int j = 0; j < 10000; ++j) {
		double* rs = new double[n + 1];
		rs[0] = r0;
		double* rand_num = getLGM(n, seed + j * 1000);
		double* z = box_muller(rand_num, n);
		for (int i = 1; i < n + 1; ++i) {
			rs[i] = rs[i - 1] + K * (r_bar - rs[i - 1])*dt + sigma * sqrt(rs[i-1])*sqrt(dt)*z[i - 1];
		}
		double P_TS = 1000 * exp(-trapzoid_method(n/2, n, rs, dt));
		double c = exp(-trapzoid_method(0, n/2, rs, dt)) * max(P_TS - strike, 0.0);
		result[j] = c;
		delete[] rand_num;
		delete[] z;
		delete[] rs;
	}
	return calc_mean(result, 10000);
}

void cir_explicit(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t) {
	double h1 = sqrt(K*K + 2 * sigma*sigma);
	double h2 = (K + h1) / 2;
	double h3 = 2 * K*r_bar / (sigma*sigma);
	double r_S = r0 * exp(-K * (T - t)) + r_bar * (1 - exp(-K * (T - t)));
	double A_TS = pow((h1*exp(h2*(S - T))) / (h2*(exp(h1*(S - T)) - 1) + h1), h3);
	double B_TS = (exp(h1*(S - T)) - 1) / (h2*(exp(h1*(S - T)) - 1) + h1);
	double P_TS = A_TS * exp(-B_TS * r_S);

	double A_tS = pow((h1*exp(h2*(S - t))) / (h2*(exp(h1*(S - t)) - 1) + h1), h3);
	double B_tS = (exp(h1*(S - t)) - 1) / (h2*(exp(h1*(S - t)) - 1) + h1);
	double P_tS = A_tS * exp(-B_tS * r0);

	double A_tT = pow((h1*exp(h2*(T - t))) / (h2*(exp(h1*(T - t)) - 1) + h1), h3);
	double B_tT = (exp(h1*(T - t)) - 1) / (h2*(exp(h1*(T - t)) - 1) + h1);
	double P_tT = A_tT * exp(-B_tT * r0);

	double theta = h1;
	double phi = 2 * theta / (sigma*sigma*(exp(theta*(T - t)) - 1));
	double psi = (K + theta) / (sigma*sigma);
	double r_star = log(A_TS / strike*1000) / B_TS;

	double x1 = 2 * r_star*(phi + psi + B_TS);
	double p1 = 4 * K*r_bar/(sigma*sigma);
	double q1 = (2 * phi*phi*r0 * exp(theta*(T - t))) / (phi + psi + B_TS);
	double x2 = 2 * r_star*(phi + psi);
	double p2 = 4 * K*r_bar / (sigma*sigma);
	double q2 = (2 * phi*phi*r0 * exp(theta*(T - t))) / (phi + psi);
	/*cout << x1 << endl;
	cout << p1 << endl;
	cout << q1 << endl;
	cout << x2 << endl;
	cout << p2 << endl;
	cout << q2 << endl;*/

	//these two values were calculated by matlab
	double ncx1 = 0.289299446954915;
	double ncx2 = 0.286418320329034;
	double call = V*P_tS*ncx1 - strike* P_tT*ncx2;
	cout << call << endl;
}

double trapzoid_method(int a, int b, double* integrand, double dt) {
	double sum = 0;
	int i = a;
	while(i < b){
		sum = sum + (integrand[i] + integrand[i + 1])*dt/2;
		++i;
	}
	return sum;
}

