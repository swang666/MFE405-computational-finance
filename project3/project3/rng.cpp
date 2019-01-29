#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "rng.h"
#include <cmath>
#include <iomanip>

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

double* question1(int64_t seed) {
	static double out[4];
	double *rand_num;
	double *rand_num2;
	//double *rand_num2 = getLGM(2000, seed + 200);
	double *z1, *z2;
	//double *z2 = box_muller(rand_num2, 2000);
	double x0 = 1;
	double y0 = 3.0 / 4.0;
	double result[1000];
	double result2[1000];
	double result3[1000];
	double result4 = 0;
	for (int i = 0; i < 1000; ++i) {
		double XT_2[2001];
		double YT_3[3001];
		XT_2[0] = x0;
		YT_3[0] = y0;
		rand_num = getLGM(2002, seed + i * 200);
		z1 = box_muller(rand_num, 2002);
		rand_num2 = getLGM(3002, seed + 1234 + i * 200);
		z2 = box_muller(rand_num2, 3002);
		double temp;
		double dt = 0.001;
		for (int j = 1; j < 2001; ++j) {
			XT_2[j] = XT_2[j - 1] + (0.2 - 0.5*XT_2[j - 1])* dt + 2.0 / 3.0 * (sqrt(j * dt)*z1[j] - sqrt((j - 1)*dt)*z1[j - 1]);
		}
		for (int j = 1; j < 3001; ++j) {
			YT_3[j] = YT_3[j - 1] + (2.0 / (1 + (j-1) * dt) * YT_3[j - 1] + (1 + pow((j-1)*dt, 3)) / 3) * dt + (1 + pow((j-1)*dt, 3)) / 3 * (sqrt(j * dt)*z2[j] - sqrt((j - 1)*dt)*z2[j - 1]);
		}

		if (YT_3[2000] > 5) {
			++result4;
		}
		result[i] = cbrt(XT_2[2000]);
		result2[i] = YT_3[3000];
		if (XT_2[2000] > 1) {
			result3[i] = XT_2[2000] * YT_3[2000];
		}
		else {
			result3[i] = 0;
		}


	}
	out[0] = result4/1000;
	out[1] = calc_mean(result, 1000);
	out[2] = calc_mean(result2, 1000);
	out[3] = calc_mean(result3, 1000);
	return out;
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
	double *z = box_muller(rand_num, 10000);
	double* w1 = wiener_process(z, 10000, 5);
	double result[16];
	double Ea1[10000];
	double control[10000];
	for (int i = 0; i < 10000; ++i) {
		Ea1[i] = w1[i] * w1[i] + sin(w1[i]);
		control[i] = w1[i] * w1[i];
	}
	result[0] = calc_mean(Ea1, 10000);
	result[8] = cov(Ea1, Ea1, 10000);
	double cov_w1_Ea1 = cov(Ea1, control, 10000);
	double var_w1 = cov(control, control, 10000);
	double gamma_1 = cov_w1_Ea1 / var_w1;

	double *Ea1_b = vector_minus_v(Ea1, 10000, vector_mult(vector_minus(control, 10000, 5), 10000, gamma_1), false);
	result[4] = calc_mean(Ea1_b, 10000);
	result[9] = cov(Ea1_b, Ea1_b, 10000);
	double* w2 = wiener_process(z, 10000, 0.5);
	double* w3 = wiener_process(z, 10000, 3.2);
	double* w4 = wiener_process(z, 10000, 6.5);
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
		con2[i] = 0.03515 * pow(w2[i], 4) - 0.6 * pow(w2[i],2) + exp(0.25);
		con3[i] = 0.1356 * pow(w3[i], 4) - 2.34 * pow(w3[i], 2) + exp(1.6);
		con4[i] = 0.706 * pow(w4[i], 4) - 12.19 * pow(w4[i], 2) + exp(3.25);
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
	result[10] = cov(Ea2, Ea2, 10000);
	result[2] = calc_mean(Ea3, 10000);
	result[11] = cov(Ea3, Ea3, 10000);
	result[3] = calc_mean(Ea4, 10000);
	result[12] = cov(Ea4, Ea4, 10000);
	double* Ea2_b = vector_minus_v(Ea2, 10000, vector_mult(vector_minus(con2, 10000, 0.03515*3*0.5*0.5 - 0.6*0.5+ exp(0.25)), 10000, gamma_2), false);
	double* Ea3_b = vector_minus_v(Ea3, 10000, vector_mult(vector_minus(con3, 10000, 0.1356 * 3 * 3.2*3.2 - 2.34*3.2 + exp(1.6)), 10000, gamma_3), false);
	double* Ea4_b = vector_minus_v(Ea4, 10000, vector_mult(vector_minus(con4, 10000, 0.706 * 3 * 6.5*6.5 - 12.19*6.5 + exp(3.25)), 10000, gamma_4), false);
	result[5] = calc_mean(Ea2_b, 10000);
	result[13] = cov(Ea2_b, Ea2_b, 10000);
	result[6] = calc_mean(Ea3_b, 10000);
	result[14] = cov(Ea3_b, Ea3_b, 10000);
	result[7] = calc_mean(Ea4_b, 10000);
	result[15] = cov(Ea4_b, Ea4_b, 10000);
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
	double out[5];
	double* rand_num = getLGM(10000, seed);
	double *z = box_muller(rand_num, 10000);
	double* p1 = euro_call(0.04, 0.2, 88, 5, 100, z, 10000, false);
	out[0] = exp(-0.04 * 5)* calc_mean(p1, 10000);
	out[3] = exp(-0.04 * 10) * cov(p1, p1, 10000);
	out[1] = black_schole(0.04, 0.2, 88, 5, 100);
	double* p2 = euro_call(0.04, 0.2, 88, 5, 100, z, 10000, true);
	out[2] = exp(-0.04 * 5)* calc_mean(p2, 10000);
	out[4] = exp(-0.04 * 10) * cov(p2, p2, 10000);
	return out;
}

double* euro_call(double r, double sigma, double S0, double T, double X, double *nums, int n, bool antithetic) {
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
	
	return result;
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

double* geometric_brownian_motion(double r, double sigma, double S0, double T, double* nums, int n) {
	double* W_T = wiener_process(nums, n, T);
	double* result = new double[n];
	for (int i = 0; i < n; ++i) {
		result[i] = S0 * exp((r - sigma * sigma / 2)*T + sigma * W_T[i]);
	}
	return result;
}

void question5(int64_t seed, double* parta, double** partb, double* partc, double** partd) {
	double out[4];
	double* rand_num = getLGM(10000, seed);
	double *z = box_muller(rand_num, 10000);
	for (int i = 0; i < 10; ++i) {
		parta[i] = calc_mean(geometric_brownian_motion(0.04, 0.18, 88, i + 1, z, 10000), 10000);
		partc[i] = calc_mean(geometric_brownian_motion(0.04, 0.35, 88, i + 1, z, 10000), 10000);
	}
	for (int i = 0; i < 6; ++i) {
		double* rand_num2 = getLGM(1000, seed+i*12345);
		double* z2 = box_muller(rand_num2, 1000);
		partb[i] = brownian_motion_path(0.04, 0.18, 88, 10, z2, 1000);
		partd[i] = brownian_motion_path(0.04, 0.35, 88, 10, z2, 1000);
	}
}

double* brownian_motion_path(double r, double sigma, double S0, double T, double* nums, int n) {
	double increment = T / (double)n;
	double* result = new double[n];
	result[0] = S0;
	for (int i = 1; i < n; ++i) {
		result[i] = result[i - 1] * exp(sigma*sqrt(increment)*nums[i - 1] + r * increment);
	}
	return result;
}

double myfunc(double x) {
	return sqrt(1 - x * x);
}

double myfunc2(double x) {
	return sqrt(1 +x);
}

double euler_approx(double t, int n, double x0) {
	double result = 0;
	double dt = t / (double)n;
	for (int i = 1; i < n; ++i) {
		result = result + myfunc(x0)*dt;
		x0 = x0 + dt;
	}
	return result;
}

double* question6(int64_t seed) {
	double out[5];
	out[0] = 4 * euler_approx(1, 10000, 0);
	double* rand_num = getLGM(10000, seed);
	double* mc = monte_carlo(rand_num, 10000, myfunc);
	out[3] = cov(mc, mc, 10000);
	out[1] = 4 * calc_mean(mc, 10000);
	double rand_num2[10000];
	for (int i = 0; i < 10000; ++i) {
		rand_num2[i] = 1-pow(rand_num[i], 2.0/3.0)*pow(3.0/2.0,2.0/3.0);
	}
	double* mc2 = monte_carlo(rand_num2, 10000, myfunc2);
	out[2] = 4 * calc_mean(mc2, 10000);
	out[4] = cov(mc2, mc2, 10000);
	return out;
}

double* monte_carlo(double* nums, int n, double (*f)(double)) {
	double* result = new double[n];
	for (int i = 0; i < n; ++i) {
		if (nums[i] < 0 || nums[i] > 1) {
			result[i] = 0;
		}
		else {
			result[i] = (*f)(nums[i]);
		}
		
	}
	return result;
}