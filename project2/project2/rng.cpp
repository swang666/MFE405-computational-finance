#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "rng.h"
#include <cmath>
#include <iomanip>

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
	double pi = atan(1) * 4;
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

double question1(int64_t seed) {
	double *rand_num = getLGM(1000, seed);
	double *rand_num2 = getLGM(1000, seed + 200);
	double *z1 = box_muller(rand_num, 1000);
	double *z2 = box_muller(rand_num2, 1000);
	double *y = generate_bivariate(z1, z2, 1000, -0.7);
	double mean_x = calc_mean(z1, 1000);
	double mean_y = calc_mean(y, 1000);
	double * x_minus_xbar = vector_minus(z1, 1000, mean_x);
	double * y_minus_ybar = vector_minus(y, 1000, mean_y);
	double numerator = 0;
	double denom1 = 0;
	double denom2 = 0;
	for (int i = 0; i < 1000; ++i) {
		numerator = numerator + x_minus_xbar[i] * y_minus_ybar[i];
		denom1 = denom1 + x_minus_xbar[i] * x_minus_xbar[i];
		denom2 = denom2 + y_minus_ybar[i] * y_minus_ybar[i];
	}
	return numerator / sqrt(denom1 * denom2);
}

double * generate_bivariate(double * z1,double * z2, int n, double rho) {
	double* y = new double[n];
	for (int i = 0; i < n; ++i) {
		y[i] = rho * z1[i] + sqrt(1 - pow(rho, 2)) * z2[i];
	}
	return y;
}

double question2(int64_t seed) {
	double *rand_num = getLGM(10000, seed);
	double *rand_num2 = getLGM(10000, seed + 200);
	double *z1 = box_muller(rand_num, 10000);
	double *z2 = box_muller(rand_num2, 10000);
	double *y = generate_bivariate(z1, z2, 10000, 0.6);
	double result[10000];
	double temp;
	for (int i = 0; i < 10000; ++i) {
		temp = pow(z1[i], 3) + sin(y[i]) + z1[i] * z1[i] * y[i];
		if (temp > 0) {
			result[i] = temp;
		}
		else {
			result[i] = 0;
		}
	}
	return calc_mean(result, 10000);
}

double * wiener_process(double* nums, int n, double t) {
	double* w = new double[n];
	for (int i = 0; i < n; ++i) {
		w[i] = sqrt(t)*nums[i];
	}
	return w;
}

double* question3(int64_t seed) {
	double* rand_num = getLGM(10000, seed);
	double* rand_num2 = getLGM(10000, seed + 234);
	double *z = box_muller(rand_num, 10000);
	double* w1 = wiener_process(z, 10000, 5);
	double result[8];
	double Ea1[10000];
	double control[10000];
	for (int i = 0; i < 10000; ++i) {
		Ea1[i] = w1[i] * w1[i] + sin(w1[i]);
		control[i] = w1[i] * w1[i];
		rand_num2[i] = 1 - rand_num[i];
	}
	result[0] = calc_mean(Ea1, 10000);
	double *z2 = box_muller(rand_num2, 10000);
	double cov_w1_Ea1 = cov(Ea1, control, 10000);
	double var_w1 = cov(control, control, 10000);
	double gamma_1 = cov_w1_Ea1 / var_w1;

	result[4] = calc_mean(vector_minus_v(Ea1, 10000, vector_mult(vector_minus(control, 10000, 5), 10000, gamma_1), false), 10000);

	double pi = atan(1) * 4;

	double* w2 = wiener_process(z, 10000, 0.5);
	double* w2_p = wiener_process(z2, 10000, 0.5);
	double* w3 = wiener_process(z, 10000, 3.2);
	double* w3_p = wiener_process(z2, 10000, 3.2);
	double* w4 = wiener_process(z, 10000, 6.5);
	double* w4_p = wiener_process(z2, 10000, 6.5);
	double Ea2[10000];
	double Ea3[10000];
	double Ea4[10000];
	double con2[10000];
	double con3[10000];
	double con4[10000];
	for (int i = 0; i < 10000; ++i) {
		Ea2[i] = exp(0.5/2) * cos(w2[i]);
		Ea3[i] = exp(3.2 / 2) * cos(w3[i]);
		Ea4[i] = exp(6.5 / 2) * cos(w4[i]);
		con2[i] = cos(w2[i]);
		con3[i] = cos(w3[i]);
		con4[i] = cos(w4[i]);
	}
	double cov_con2_Ea2 = cov(Ea2, con2, 10000);
	double var_con2 = cov(con2, con2, 10000);
	double gamma_2 = cov_con2_Ea2 / var_con2;
	double cov_con3_Ea3 = cov(Ea3, con3, 10000);
	double var_con3 = cov(con3, con3, 10000);
	double gamma_3 = cov_con3_Ea3 / var_con3;
	double cov_con4_Ea4 = cov(Ea4, con4, 10000);
	double var_con4 = cov(con4, con4, 10000);
	double gamma_4 = cov_con4_Ea4 / var_con4;

	result[1] = calc_mean(Ea2, 10000);
	result[2] = calc_mean(Ea3, 10000);
	result[3] = calc_mean(Ea4, 10000);
	result[5] = calc_mean(vector_minus_v(Ea2, 10000, vector_mult(vector_minus(con2, 10000, exp(-0.5 / 2)), 10000, gamma_2), false), 10000);
	result[6] = calc_mean(vector_minus_v(Ea3, 10000, vector_mult(vector_minus(con3, 10000, exp(-3.2 / 2)), 10000, gamma_3), false), 10000);
	result[7] = calc_mean(vector_minus_v(Ea4, 10000, vector_mult(vector_minus(con4, 10000, exp(-6.5 / 2)), 10000, gamma_4), false), 10000);
	return result;
}

double cov(double* X, double* Y, int n) {
	double mean_x = calc_mean(X, n);
	double* x_minus_xbar = vector_minus(X, n, mean_x);
	double mean_y = calc_mean(Y, n);
	double* y_minus_ybar = vector_minus(Y, n, mean_y);
	double result = 0;
	for (int i = 0; i < n; ++i) {
		result = result + x_minus_xbar[i] * y_minus_ybar[i];
	}
	result = result / (n - 1);
	return result;
}

double* question4(int64_t seed) {
	//a
	double out[4];
	double* rand_num = getLGM(10000, seed);
	double *z = box_muller(rand_num, 10000);
	out[0] = euro_call(0.04, 0.2, 88, 5, 100, z, 10000, false);
	out[1] = black_schole(0.04, 0.2, 88, 5, 100);
	out[2] = euro_call(0.04, 0.2, 88, 5, 100, z, 10000, true);
	return out;
}

double euro_call(double r, double sigma, double S0, double T, double X, double *nums, int n, bool antithetic) {
	double* W_T = wiener_process(nums, n, T);
	double* result = new double[n];
	double temp;
	if (antithetic == false) {
		for (int i = 0; i < n; ++i) {
			temp = S0 * exp((r - sigma * sigma / 2)*T + sigma * W_T[i]) - X;
			if (temp > 0) {
				result[i] = temp;
			}
			else {
				result[i] = 0;
			}
		}
	}
	else {
		double* z2 = new double[n];
		double temp2;
		for (int i = 0; i < n; ++i) {
			z2[i] = -nums[i];
		}
		double* W_T2 = wiener_process(z2, n, T);
		for (int i = 0; i < n; ++i) {
			temp = S0 * exp((r - sigma * sigma / 2)*T + sigma * W_T[i]) - X;
			temp2 = S0 * exp((r - sigma * sigma / 2)*T + sigma * W_T2[i]) - X;
			if (temp > 0 && temp2 > 0) {
				result[i] = (temp + temp2)/2;
			}
			else if (temp > 0){
				result[i] = temp/2;
			}
			else if (temp2 > 0) {
				result[i] = temp2 / 2;
			}
			else {
				result[i] = 0;
			}
		}
	}
	
	return exp(-r * T)*calc_mean(result, n);
}

double normalCDF(double x) 
{
	return std::erfc(-x / std::sqrt(2)) / 2;
}

double black_schole(double r, double sigma, double S0, double T, double X) {
	double d1 = (log(S0 / X) + (r + sigma * sigma / 2)*T) / (sigma*sqrt(T));
	double d2 = d1 - sigma * sqrt(T);
	return S0 * normalCDF(d1) - X * exp(-r * T) * normalCDF(d2);
}