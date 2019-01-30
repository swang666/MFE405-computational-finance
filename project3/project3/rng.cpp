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
		rand_num = getLGM(2000, seed + i * 200);
		z1 = box_muller(rand_num, 2000);
		rand_num2 = getLGM(3000, seed + 1234 + i * 200);
		z2 = box_muller(rand_num2, 3000);
		double dt = 0.001;
		for (int j = 1; j < 2001; ++j) {
			XT_2[j] = XT_2[j - 1] + (0.2 - 0.5*XT_2[j - 1])* dt + 2.0 / 3.0 * sqrt(dt)*z1[j - 1];
		}
		for (int j = 1; j < 3001; ++j) {
			YT_3[j] = YT_3[j - 1] + (2.0 / (1 + (j-1) * dt) * YT_3[j - 1] + (1 + pow((j-1)*dt, 3)) / 3) * dt + (1 + pow((j-1)*dt, 3)) / 3 * sqrt(dt)*z2[j - 1];
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


double* question2(int64_t seed) {
	static double out[2];
	double *rand_num;
	double *rand_num2;
	//double *rand_num2 = getLGM(2000, seed + 200);
	double *w, *z;
	//double *z2 = box_muller(rand_num2, 2000);
	double x0 = 1;
	double y0 = 1;
	double result[1000];
	double result2[1000];
	for (int i = 0; i < 1000; ++i) {
		double XT_3[3001];
		double log_YT_3[3001];
		XT_3[0] = x0;
		log_YT_3[0] = log(y0);
		rand_num = getLGM(3000, seed + i * 200);
		w = box_muller(rand_num, 3000);
		rand_num2 = getLGM(3000, seed + 12345 + i * 200);
		z = box_muller(rand_num2, 3000);
		double dt = 0.001;
		for (int j = 1; j < 3001; ++j) {
			XT_3[j] = XT_3[j-1] + 0.25 * XT_3[j - 1] * dt + 1.0 / 3.0 * XT_3[j - 1] * sqrt(dt) * w[j - 1] - 0.75 * XT_3[j - 1] * sqrt(dt) * z[j - 1] + 0.5 / 3.0 * XT_3[j - 1] / 3.0 * (sqrt(dt) * w[j - 1] * sqrt(dt) * w[j - 1] - dt) - 0.5 * 0.75 * XT_3[j - 1] * 0.75 * (sqrt(dt) * z[j - 1] * sqrt(dt) * z[j - 1] - dt);
			log_YT_3[j] = log_YT_3[j - 1] - 0.08  * dt + 1.0 / 3.0 * sqrt(dt) * w[j - 1] + 0.75 * sqrt(dt) * z[j - 1];
		}
		result[i] = cbrt(1 + XT_3[3000]);
		result2[i] = cbrt(1 + exp(log_YT_3[3000]));
		
	}
	out[0] = calc_mean(result, 1000);
	out[1] = calc_mean(result2, 1000);
	return out;
}

double * wiener_process(double* nums, int n, double t) {
	double* w = new double[n];
	for (int i = 0; i < n; ++i) {
		w[i] = sqrt(t)*nums[i];
	}
	return w;
}

void question3(int64_t seed, double S0, double T, double X, double r, double sigma, double* prices, double* prices_2, double* parta, double** partc, double** partc_2) {
	double p1 = monte_carlo_euro_call(S0, T, X, r, sigma, seed);
	double p2 = black_schole(r, sigma, S0, T, X);
	parta[0] = p1;
	parta[1] = p2;
	for (int i = 0; i < 11; ++i) {
		partc[i] = greeks(S0 + i, T, r, X, sigma);
		partc_2[i] = greeks_2(S0 + i, T, r, X, sigma, seed);
	}
	double d = 0.01;
	for (int i = 0; i < 1100; ++i) {
		prices[i] = monte_carlo_euro_call(S0 + i * d, T, X, r, sigma, seed);
		prices_2[i] = black_schole(r, sigma, S0 + i*d, T, X);
	}
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
	double* out = new double[3];
	double reflec = heston_model(-0.6, 0.03, 48, 0.05, 0.42, 5.8, 0.0625, 0.5, 50,0, seed);
	double p_trun = heston_model(-0.6, 0.03, 48, 0.05, 0.42, 5.8, 0.0625, 0.5, 50,1, seed);
	double f_trun = heston_model(-0.6, 0.03, 48, 0.05, 0.42, 5.8, 0.0625, 0.5, 50,2, seed);
	out[0] = exp(-0.03 * 0.5) * reflec;
	out[1] = exp(-0.03 * 0.5) * p_trun;
	out[2] = exp(-0.03 * 0.5) * f_trun;
	return out;
}

double monte_carlo_euro_call(double S0, double T, double X, double r, double sigma, int64_t seed) {
	double result[1000];
	double* rand_num = getLGM(1000, seed);
	double* z = box_muller(rand_num, 1000);
	double* W_T = wiener_process(z, 1000, T);
	double W_T2[1000];
	for (int i = 0; i < 1000; ++i) {
		W_T2[i] = -W_T[i];
	}
	double temp, temp2;
	for (int i = 0; i < 1000; ++i) {
		temp = S0 * exp((r - sigma * sigma / 2)*T + sigma * W_T[i]) - X;
		temp2 = S0 * exp((r - sigma * sigma / 2)*T + sigma * W_T2[i]) - X;
		if (temp > 0 && temp2 > 0) {
			result[i] = (temp + temp2) / 2;
		}
		else if (temp > 0) {
			result[i] = temp / 2;
		}
		else if (temp2 > 0) {
			result[i] = temp2 / 2;
		}
		else {
			result[i] = 0;
		}
	}
	return exp(-r * T)* calc_mean(result, 1000);
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

double pdf(double x) {
	return 1 / sqrt(2 * pi) * exp(-x * x / 2);
}

double sim_normalCDF(double x) {
	double d1 = 0.049867347, d2 = 0.0211410061, d3 = 0.0032776263, d4 = 0.0000380036, d5 = 0.0000488906, d6 = 0.000005383;
	if (x >= 0) {
		return (1 - 0.5*pow((1 + d1 * x + d2 * pow(x, 2) + d3 * pow(x, 3) + d4 * pow(x, 4) + d5 * pow(x, 5) + d6 * pow(x, 6)),-16));
	}
	else {
		return (0.5*pow((1 + d1 * -x + d2 * pow(-x, 2) + d3 * pow(-x, 3) + d4 * pow(-x, 4) + d5 * pow(-x, 5) + d6 * pow(-x, 6)), -16));
	}
}

double normalCDF(double x) 
{
	return std::erfc(-x / std::sqrt(2)) / 2;
}

double black_schole(double r, double sigma, double S0, double T, double X) {
	double d1 = (log(S0 / X) + (r + sigma * sigma / 2)*T) / (sigma*sqrt(T));
	double d2 = d1 - sigma * sqrt(T);
	return S0 * sim_normalCDF(d1) - X * exp(-r * T) * sim_normalCDF(d2);
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

double* greeks(double S0, double T, double r, double X, double sigma) {
	double* out = new double[5];
	double d = 0.01;
	double delta = (black_schole(r, sigma, S0+d, T, X) - black_schole(r, sigma, S0, T, X))/d;
	out[0] = delta;
	double gamma = (black_schole(r, sigma, S0 + d, T, X) - 2 * black_schole(r, sigma, S0, T, X) + black_schole(r, sigma, S0 - d, T, X)) / (d*d);
	out[1] = gamma;
	double theta = (black_schole(r, sigma, S0, T +d, X) - black_schole(r, sigma, S0, T, X)) / d;
	out[2] = theta;
	double vega = (black_schole(r, sigma+ d*0.01, S0, T, X) - black_schole(r, sigma, S0, T, X)) / (d*0.01);
	out[3] = vega;
	double rho = (black_schole(r+ d*0.01, sigma, S0, T, X) - black_schole(r, sigma, S0, T, X)) / (d*0.01);
	out[4] = rho;
	return out;

}

double* greeks_2(double S0, double T, double r, double X, double sigma, int64_t seed) {
	double* out = new double[5];
	double d = 0.01;
	double delta = (monte_carlo_euro_call(S0+d, T, X, r, sigma, seed) - monte_carlo_euro_call(S0, T, X, r, sigma, seed)) / d;
	out[0] = delta;
	double gamma = (monte_carlo_euro_call(S0+d, T, X, r, sigma, seed) - 2 * monte_carlo_euro_call(S0, T, X, r, sigma, seed) + monte_carlo_euro_call(S0-d, T, X, r, sigma, seed)) / (d*d);
	out[1] = gamma;
	double theta = (monte_carlo_euro_call(S0, T+d, X, r, sigma, seed) - monte_carlo_euro_call(S0, T, X, r, sigma, seed)) / d;
	out[2] = theta;
	double vega = (monte_carlo_euro_call(S0, T, X, r, sigma+ d*0.01, seed) - monte_carlo_euro_call(S0, T, X, r, sigma, seed)) / (d*0.01);
	out[3] = vega;
	double rho = (monte_carlo_euro_call(S0, T, X, r + d*0.01, sigma, seed) - monte_carlo_euro_call(S0, T, X, r, sigma, seed)) / (d*0.01);
	out[4] = rho;
	return out;
}

double heston_model(double rho, double r, double S0, double V0, double sigma, double alpha, double beta, double T, double K, int type, int64_t seed) {
	double out;
	double* rand_num;
	double* rand_num2;
	double* w1;
	double* w2;
	double dt = 0.01;
	int size = T / dt;
	double* result = new double[1000];
	double* V = new double[size+1];
	V[0] = V0;
	double* S = new double[size + 1];
	S[0] = S0;
	switch (type)
	{
	case 0:
		for (int i = 0; i < 1000; ++i) {
			rand_num = getLGM(size * 2, seed + i * 200);
			rand_num2 = getLGM(size * 2, seed + 1234 + i * 200);
			w1 = box_muller(rand_num, size * 2);
			double* t = box_muller(rand_num2, size * 2);
			w2 = generate_bivariate(w1, t, size * 2, rho);
			for (int j = 1; j < size+1; ++j) {
				V[j] = abs(V[j - 1]) + alpha * (beta - abs(V[j - 1])) * dt + sigma * sqrt(abs(V[j - 1])) * sqrt(dt)*w1[j];
				S[j] = S[j - 1] + r * S[j - 1] * dt + sqrt(abs(V[j - 1]))*S[j - 1] * sqrt(dt)* w2[j];
			}
			result[i] = max(S[size] - K, 0.0);
		}
		out = calc_mean(result, 1000);
		break;
	case 1:
		for (int i = 0; i < 1000; ++i) {
			rand_num = getLGM(size * 2, seed + i * 200);
			rand_num2 = getLGM(size * 2, seed + 1234 + i * 200);
			w1 = box_muller(rand_num, size * 2);
			double* t = box_muller(rand_num2, size * 2);
			w2 = generate_bivariate(w1, t, size * 2, rho);
			for (int j = 1; j < size + 1; ++j) {
				V[j] = V[j - 1] + alpha * (beta - V[j - 1]) * dt + sigma * sqrt(max(V[j-1],0.0)) * sqrt(dt)*w1[j];
				S[j] = S[j - 1] + r * S[j - 1] * dt + sqrt(max(V[j - 1], 0.0))*S[j - 1] * sqrt(dt)*w2[j];
			}
			result[i] = max(S[size] - K, 0.0);
		}
		out = calc_mean(result, 1000);
		break;
	default:
		for (int i = 0; i < 1000; ++i) {
			rand_num = getLGM(size * 2, seed + i * 200);
			rand_num2 = getLGM(size * 2, seed + 1234 + i * 200);
			w1 = box_muller(rand_num, size * 2);
			double* t = box_muller(rand_num2, size * 2);
			w2 = generate_bivariate(w1, t, size * 2, rho);
			for (int j = 1; j < size + 1; ++j) {
				V[j] = V[j - 1] + alpha * (beta - max(V[j - 1], 0.0)) * dt + sigma * sqrt(max(V[j - 1], 0.0)) * sqrt(dt)*w1[j];
				S[j] = S[j - 1] + r * S[j - 1] * dt + sqrt(max(V[j - 1], 0.0))*S[j - 1] * sqrt(dt)*w2[j];
			}
			result[i] = max( S[size] - K , 0.0);
		}
		out = calc_mean(result, 1000);
		break;
	}
	return out;
}