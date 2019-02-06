#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "binom.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;

double** question1() {
	double T = 0.5;
	double dt, c, u, d, p;
	double r = 0.05;
	double sigma = 0.24;
	double S0 = 32;
	double K = 30;
	int N[7] = { 10, 20, 40, 80, 100, 200, 500 };
	double** out = new double*[4];
	for (int j = 0; j < 4; ++j) {
		out[j] = new double[7];
	}
	for (int i = 0; i < 7; ++i) {
			dt = T / N[i];
			//cout << N[i] << endl;
			c = 0.5*(exp(-r * dt) + exp((r + sigma * sigma)*dt));
			d = c - sqrt(c*c - 1);
			//cout << d << endl;
			u = 1 / d;
			//cout << u << endl;
			p = (exp(r*dt) - d) / (u - d);
			//cout << p << endl;
			out[0][i] = binom_method2(S0, K, u, d, p, N[i], T, r,0,0);
	}
	for (int i = 0; i < 7; ++i) {
		dt = T / N[i];
		//cout << N[i] << endl;

		d = exp(r*dt)*(1 - sqrt(exp(sigma*sigma*dt) - 1));
		//cout << d << endl;
		u = exp(r*dt)*(1 + sqrt(exp(sigma*sigma*dt) - 1));

		//cout << u << endl;
		p = 0.5;
		//cout << p << endl;
		out[1][i] = binom_method2(S0, K, u, d, p, N[i], T, r,0,0);
	}
	for (int i = 0; i < 7; ++i) {
		dt = T / N[i];
		//cout << N[i] << endl;

		d = exp((r - sigma * sigma / 2)*dt - sigma * sqrt(dt));
		//cout << d << endl;
		u = exp((r - sigma * sigma / 2)*dt + sigma * sqrt(dt));

		//cout << u << endl;
		p = 0.5;
		//cout << p << endl;
		out[2][i] = binom_method2(S0, K, u, d, p, N[i], T, r,0,0);
	}
	for (int i = 0; i < 7; ++i) {
		dt = T / N[i];
		//cout << N[i] << endl;

		d = exp(-sigma * sqrt(dt));
		//cout << d << endl;
		u = exp(sigma*sqrt(dt));

		//cout << u << endl;
		p = 0.5 + 0.5*((r-sigma*sigma/2)*sqrt(dt)/sigma);
		//cout << p << endl;
		out[3][i] = binom_method2(S0, K, u, d, p, N[i], T, r,0,0);
	}
	return out;
}

double question2() {
	double r = 0.02;
	double sigma = 0.202453;
	double K = 1265;
	double S0 = 1145.99;
	double T = 1;
	int n = 1000;
	double dt = T / n;
	double d, u, p;
	d = exp((r - sigma * sigma / 2)*dt - sigma * sqrt(dt));
	//cout << d << endl;
	u = exp((r - sigma * sigma / 2)*dt + sigma * sqrt(dt));

	//cout << u << endl;
	p = 0.5;
	return binom_method2(S0, K, u, d, p, n, T, r,0,0);
}

double** question3() {
	double** out = new double*[6];
	double r = 0.03;
	double sigma = 0.2;
	double K = 50;
	double S = 49;
	double S0 = 20;
	double T = 0.3846;
	double T0 = 0.01;
	int n1 = (80 - 20) / 2;
	double d = 0.01;
	double dt1 = T / n1;
	double dt2 = (T + 0.01) / n1;
	double* udp4 = JR_Btree(r, sigma, dt2);
	double* udp1 = JR_Btree(r, sigma, dt1);
	double* udp2 = JR_Btree(r, sigma + d * 0.01, dt1);
	double* udp3 = JR_Btree(r+ d* 0.01, sigma, dt1);
	for (int i = 0; i < 5; ++i) {
		out[i] = new double[n1 + 1];
	}
	for (int i = 0; i < n1+1; ++i) {
		//delta
		out[0][i] = (binom_method2(S0+(i+1)*2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0) - binom_method2(S0 + i * 2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0))/2;
		//theta
		out[1][i] = (binom_method2(S0 + i * 2, K, udp4[0], udp4[1], udp4[2], n1, T+0.01, r,0,0) - binom_method2(S0 + i * 2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0)) / d;
		//gamma
		out[2][i] = (binom_method2(S0 + (i+1) * 2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0) - 2* binom_method2(S0 + i * 2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0) + binom_method2(S0 + (i-1) * 2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0)) / (4);
		//vega
		out[3][i] = (binom_method2(S0 + i * 2, K, udp2[0], udp2[1], udp2[2], n1, T, r,0,0) - binom_method2(S0 + i * 2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0)) / (d * 0.01);
		//rho
		out[4][i] = (binom_method2(S0 + i * 2, K, udp3[0], udp3[1], udp3[2], n1, T, r+ d*0.01,0,0) - binom_method2(S0 + i * 2, K, udp1[0], udp1[1], udp1[2], n1, T, r,0,0)) / (d * 0.01);
	}
	int n2 = T / 0.01;
	out[5] = new double[n2];
	for (int i = 0; i < n2; ++i) {
		double t_dt1 = (T0 + i * 0.01) / n2;
		double* t_udp1 = JR_Btree(r, sigma, t_dt1);
		out[5][i] = (binom_method2(S+0.01, K, t_udp1[0], t_udp1[1], t_udp1[2], n2, T0 + i * 0.01, r ,0,0) - binom_method2(S, K, t_udp1[0], t_udp1[1], t_udp1[2], n2, T0 + i*0.01, r,0,0)) / (0.01);
	}
	return out;
}

double** question4() {
	double** out = new double*[2];
	out[0] = new double[11];
	out[1] = new double[11];
	double S0 = 80;
	double K = 100;
	double r = 0.05;
	double sigma = 0.3;
	double T = 1;
	int n = 100;
	double* udp = JR_Btree(r, sigma, T / n);
	for (int i = 0; i < 11; ++i) {
		out[0][i] = binom_method2(S0+ i*4, K, udp[0], udp[1], udp[2], n, T, r, 0, 1);
		out[1][i] = binom_method2(S0+ i*4, K, udp[0], udp[1], udp[2], n, T, r, 1, 1);
	}
	return out;
}

double binom_method2(double S0, double K, double u, double d, double p, int n, double T, double r, int type, int put) {
	double h = T / n;
	double** prices = new double*[n+1];
	double** payoffs = new double*[n + 1];
	double** p_no = new double*[n];
	double** p_ex = new double*[n];
	for (int i = n; i >= 0; --i) {
		prices[i] = new double[n + 1];
		payoffs[i] = new double[n + 1];

		for (int j = 0; j <= i; ++j) {
			prices[i][j] = S0 * pow(u, (i - j))*pow(d , j);
			if (i == n) {
				//initialize final payoff
				if (put == 0)
					payoffs[i][j] = max(prices[i][j] - K, 0.0);
				else
					payoffs[i][j] = max(-prices[i][j] + K, 0.0);
			}
			else {
				if (type == 0) {
					payoffs[i][j] = exp(-r * h)*(payoffs[i + 1][j] * p + payoffs[i + 1][j + 1] * (1-p));
					//payoffs(j, i) = exp(-r * h)* (payoffs(j, i + 1)*(exp(r*h) - d) / (u - d) + payoffs(j + 1, i + 1)*(u - exp(r*h)) / (u - d));
				}
				else {
					p_no[i] = new double[n];
					p_ex[i] = new double[n];
					p_no[i][j] = exp(-r * h)*(payoffs[i + 1][j] * p + payoffs[i + 1][j + 1] * (1 - p));
					if (put == 0)
						p_ex[i][j] = max(prices[i][j] - K, 0.0);
					else
						p_ex[i][j] = max(-prices[i][j] + K, 0.0);
					payoffs[i][j] = max(p_no[i][j], p_ex[i][j]);
				}
				
			}
			//cout << i << " " << payoffs[i][j] << endl;
		}
	}
	return payoffs[0][0];
}

double* JR_Btree(double r, double sigma, double dt) {
	double* out = new double[3];
	out[0] = exp((r - sigma * sigma / 2)*dt + sigma * sqrt(dt));
	out[1] = exp((r - sigma * sigma / 2)*dt - sigma * sqrt(dt));
	out[2] = 0.5;
	return out;
}