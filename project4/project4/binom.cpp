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
			out[0][i] = binom_method2(S0, K, u, d, p, N[i], T, r);
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
		out[1][i] = binom_method2(S0, K, u, d, p, N[i], T, r);
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
		out[2][i] = binom_method2(S0, K, u, d, p, N[i], T, r);
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
		out[3][i] = binom_method2(S0, K, u, d, p, N[i], T, r);
	}
	return out;
}


double binom_method2(double S0, double K, double u, double d, double p, int n, double T, double r) {
	double h = T / n;
	double** prices = new double*[n+1];
	double** payoffs = new double*[n + 1];
	for (int i = n; i >= 0; --i) {
		prices[i] = new double[n + 1];
		payoffs[i] = new double[n + 1];
		for (int j = 0; j <= i; ++j) {
			prices[i][j] = S0 * pow(u, (i - j))*pow(d , j);
			if (i == n) {
				//initialize final payoff
				payoffs[i][j] = max(prices[i][j] - K, 0.0);
				
			}
			else {
				payoffs[i][j] = exp(-r * h)*(payoffs[i + 1][j] * p + payoffs[i + 1][j + 1] * (1-p));
				//payoffs(j, i) = exp(-r * h)* (payoffs(j, i + 1)*(exp(r*h) - d) / (u - d) + payoffs(j + 1, i + 1)*(u - exp(r*h)) / (u - d));
			}
			//cout << i << " " << payoffs[i][j] << endl;
		}
	}
	return payoffs[0][0];
}
