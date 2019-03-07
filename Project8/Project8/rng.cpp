#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "rng.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>

using namespace std;

const double pi = atan(1) * 4;

default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);

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
	delete[] x;
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
	double out = calc_mean(result, 10000);
	delete[] result;
	return out;
}

double vasicek_model_semicoupon_bond(double r0, double sigma, double C, double K, double r_bar, double V, double T, double t, int64_t seed) {
	int period = T * 2;
	double* result = new double[10000];
	for (int j = 0; j < 10000; ++j) {
		double sum = 0;
		int n = 252 * T;
		double dt = (T - t) / n;
		double* rs = new double[n + 1];
		rs[0] = r0;
		double* rand_num = getLGM(n, seed + j * 100);
		double* z = box_muller(rand_num, n);
		for (int i = 1; i < n + 1; ++i) {
			rs[i] = rs[i - 1] + K * (r_bar - rs[i - 1])*dt + sigma * sqrt(dt)*z[i - 1];
		}
		for (int i = 0; i < period - 1; ++i) {
			sum = sum + C * exp(-trapzoid_method(0, n/8*(i+1), rs, dt));
		}
		sum = sum + (V+C) * exp(-trapzoid_method(0, n, rs, dt));
		result[j] = sum;
		delete[] rand_num;
		delete[] z;
		delete[] rs;
	}
	double out = calc_mean(result, 10000);
	delete[] result;
	return out;
}

double vasicek_model_zero_bond_call(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t, int64_t seed) {

	int num1 = 1000;
	int num2 = 1000;
	int n = 126;
	double dt = (S - t) / n;

	double* callprice = new double[num1];
	for (int k = 0; k < num1; ++k) {
		double* z1;
		double* Payoff = new double[num2];
		double* rT = new double[n / 2 + 1];
		rT[0] = r0;
		double* rand_num1 = getLGM(n / 2, seed*(k + 1) + 30);
		z1 = box_muller(rand_num1, n / 2);
		for (int i = 1; i < n / 2 + 1; ++i) {
			rT[i] = rT[i - 1] + K * (r_bar - rT[i - 1])*dt + sigma*sqrt(dt)*z1[i - 1];

		}
		delete[] rand_num1;
		//delete[] z1;
		for (int j = 0; j < num2; ++j) {
			
			double* rS = new double[n / 2 + 1];
			rS[0] = rT[n / 2];
			double* rand_num2 = getLGM(n / 2, (seed + k) *(j + 1) + 13);
			double* z2 = box_muller(rand_num2, n / 2);
			for (int i = 1; i < n / 2 + 1; ++i) {
				rS[i] = rS[i - 1] + K * (r_bar - rS[i - 1])*dt + sigma*sqrt(dt)*z2[i - 1];
			}
			double P_TS = V * exp(-trapzoid_method(0, n / 2, rS, dt));
			Payoff[j] = P_TS;
			delete[] rS;
			delete[] rand_num2;
			//delete[] z2;
		}

		double payoff = max(calc_mean(Payoff, num2) - strike, 0.0);
		callprice[k] = exp(-trapzoid_method(0, n / 2, rT, dt)) * payoff;
		delete[] Payoff;
		delete[] rT;
	}
	double out = calc_mean(callprice, num1);
	delete[] callprice;
	return out;
}

double cir_model_zero_bond(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t, int64_t seed) {

	int num1 = 1000;
	int num2 = 1000;
	int n = 252;
	double dt = (S - t) / n;

	double* callprice = new double[num1];
	for (int k = 0; k < num1; ++k) {
		double* Payoff = new double[num2];
		double* rT = new double[n/2 + 1];
		rT[0] = r0;
		double* rand_num1 = getLGM(n/2, seed*(k+1)+30);
		double* z1 = box_muller(rand_num1, n/2);
		for (int i = 1; i < n / 2 + 1; ++i) {
			rT[i] = rT[i - 1] + K * (r_bar - rT[i - 1])*dt + sigma * sqrt(abs(rT[i - 1]))*sqrt(dt)*z1[i-1];
		}
		delete[] rand_num1;
		delete[] z1;
		for (int j = 0; j < num2; ++j) {
			double* rS = new double[n/2 + 1];
			rS[0] = rT[n/2];
			double* rand_num2 = getLGM(n/2, (seed+k) *(j+1) + 13);
			double* z2 = box_muller(rand_num2, n/2);
			for (int i = 1; i < n / 2 + 1; ++i) {
				rS[i] = rS[i - 1] + K * (r_bar - rS[i - 1])*dt + sigma * sqrt(abs(rS[i - 1]))*sqrt(dt)*z2[i-1];
			}
			double P_TS = 1000 * exp(-trapzoid_method(0, n/2, rS, dt));
			Payoff[j] = P_TS;
			delete[] rS;
			delete[] rand_num2;
			delete[] z2;
		}

		double payoff = max(calc_mean(Payoff, num2) - strike, 0.0);
		callprice[k] = exp(-trapzoid_method(0, n / 2, rT, dt)) * payoff;

		delete[] Payoff;
		delete[] rT;
	}
	return calc_mean(callprice, num1);
}

double cir_explicit(double r0, double sigma, double K, double strike, double r_bar, double V, double T, double S, double t) {
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
	return call;
}

double G2pp(double r0, double x0, double y0, double phi0, double rho, double a, double b, double sigma, double eta, double phi_t, double strike, double V, double T, double S, double t) {
	double V_tS = sigma * sigma / (a*a)*(S - t + 2 / a * exp(-a * (S - t)) - 0.5 / a * exp(-2 * a*(S - t)) - 1.5 / a) + eta * eta / (b*b)*(S - t + 2 / b * exp(-b * (S - t)) - 0.5 / b * exp(-2 * b*(S - t)) - 1.5 / b) + 2 * rho*sigma*eta / (a*b)*(S - t + (exp(-a * (S - t)) - 1) / a + (exp(-b * (S - t)) - 1) / b - (exp(-(a + b) * (S - t)) - 1) / (a + b));
	double P_tS = exp(-(S - t)*phi_t - 0.5*V_tS);

	double V_tT = sigma * sigma / (a*a)*(T - t + 2 / a * exp(-a * (T - t)) - 0.5 / a * exp(-2 * a*(T - t)) - 1.5 / a) + eta * eta / (b*b)*(T - t + 2 / b * exp(-b * (T - t)) - 0.5 / b * exp(-2 * b*(T - t)) - 1.5 / b) + 2 * rho*sigma*eta / (a*b)*(T - t + (exp(-a * (T - t)) - 1) / a + (exp(-b * (T - t)) - 1) / b - (exp(-(a + b) * (T - t)) - 1) / (a + b));
	double P_tT = exp(-(T - t)*phi_t - 0.5*V_tT);

	double Sigma_square = sigma * sigma / (2 * a*a*a)*(1 - exp(-a * (S - T)))*(1 - exp(-a * (S - T)))*(1 - exp(-2 * a * (T - t))) + eta * eta / (2 * b*b*b)*(1 - exp(-b * (S - T)))*(1 - exp(-b * (S - T)))*(1 - exp(-2 * b * (T - t))) + 2 * rho*sigma*eta / (a*b*(a + b))*(1 - exp(-a * (S - T)))*(1 - exp(-b * (S - T)))*(1 - exp(-(a + b) * (T - t)));
	double Sigma = sqrt(Sigma_square);

	double put = -V * P_tS * normalCDF(log(strike / 1000 * P_tT / P_tS) / Sigma - 0.5*Sigma) + P_tT * strike * normalCDF(log(strike / 1000 * P_tT / P_tS) / Sigma + 0.5*Sigma);
	return put;
}

double G2pp_mc(double r0, double x0, double y0, double phi0, double rho, double a, double b, double sigma, double eta, double phi_t, double strike, double V, double T, double S, double t, int64_t seed) {
	int num1 = 1000;
	int num2 = 1000;
	int n = 252;
	double dt = (S - t) / n;
	double* putprice = new double[num1];
	for (int k = 0; k < num1; ++k) {
		double* Payoff = new double[num2];
		double* xT = new double[n/2+1];
		double* yT = new double[n/2+1];
		double* rT = new double[n / 2 + 1];
		xT[0] = x0;
		yT[0] = y0;
		rT[0] = r0;
		double* rand_num1 = getLGM(n / 2, seed *(k + 1) + 7);
		double* rand_num2 = getLGM(n / 2, seed *(k + 1) + 207);
		double* z1 = box_muller(rand_num1, n / 2);
		double* t = box_muller(rand_num2, n / 2);
		double* z2 = generate_bivariate(z1, t, n / 2, rho);
		for (int i = 1; i < n / 2 + 1; ++i) {
			xT[i] = xT[i - 1] - a * xT[i - 1] * dt + sigma * sqrt(dt)*z1[i - 1];
			yT[i] = yT[i - 1] - b * yT[i - 1] * dt + eta * sqrt(dt)*z2[i - 1];
			rT[i] = xT[i] + yT[i] + phi_t;
		}
		delete[] rand_num1;
		delete[] rand_num2;
		delete[] z1;
		delete[] t;
		delete[] z2;
		for (int j = 0; j < num2; ++j) {
			double* xS = new double[n / 2 + 1];
			double* yS = new double[n / 2 + 1];
			double* rS = new double[n / 2 + 1];
			xS[0] = xT[n / 2];
			yS[0] = yT[n / 2];
			rS[0] = rT[n / 2];
			double* rand_num1 = getLGM(n / 2, (seed+k) *(j + 1) + 7);
			double* rand_num2 = getLGM(n / 2, (seed+k) *(j + 1) + 207);
			double* z1 = box_muller(rand_num1, n / 2);
			double* t = box_muller(rand_num2, n / 2);
			double* z2 = generate_bivariate(z1, t, n / 2, rho);
			for (int i = 1; i < n / 2 + 1; ++i) {
				xS[i] = xS[i - 1] - a * xS[i - 1] * dt + sigma * sqrt(dt)*z1[i - 1];
				yS[i] = yS[i - 1] - b * yS[i - 1] * dt + eta * sqrt(dt)*z2[i - 1];
				rS[i] = xS[i] + yS[i] + phi_t;
			}
			double P_TS = V * exp(-trapzoid_method(0, n / 2, rS, dt));
			Payoff[j] = P_TS;
			delete[] rand_num1;
			delete[] rand_num2;
			delete[] z1;
			delete[] t;
			delete[] z2;
			delete[] xS;
			delete[] yS;
			delete[] rS;
		}
		double payoff = max(-calc_mean(Payoff, num2) + strike, 0.0);
		putprice[k] = exp(-trapzoid_method(0, n / 2, rT, dt)) * payoff;

		delete[] Payoff;
		delete[] rT;
		delete[] xT;
		delete[] yT;
	}
	return calc_mean(putprice, num1);
}

double * generate_bivariate(double * z1, double * z2, int n, double rho) {
	double* y = new double[n];
	for (int i = 0; i < n; ++i) {
		y[i] = rho * z1[i] + sqrt(1 - rho*rho) * z2[i];
	}
	return y;
}

double normalCDF(double x)
{
	return std::erfc(-x / std::sqrt(2)) / 2;
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

