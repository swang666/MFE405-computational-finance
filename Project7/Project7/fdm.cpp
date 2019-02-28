#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "fdm.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>
#include <Eigen>

using namespace std;
using namespace Eigen;

double EFDM(double S0, int k, int euro, int call) {
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
	int total_num_path = 2 * num_path + 1;
	VectorXd XN = VectorXd::Zero(total_num_path);
	VectorXd FN(total_num_path);
	if (call == 0) {
		// put option
		for (int i = num_path; i < 2 * num_path + 1; ++i) {
			XN(i) = log(S0) - (i - num_path) * dx;
			FN(i) = max(K - exp(XN(i)), 0.0);
		}
		for (int i = num_path - 1; i >= 0; --i) {
			XN(i) = log(S0) + (num_path - i) * dx;
			FN(i) = max(K - exp(XN(i)), 0.0);
		}
	}
	else {
		// call option
		for (int i = num_path; i < 2 * num_path + 1; ++i) {
			XN(i) = log(S0) - (i - num_path) * dx;
			FN(i) = max(-K + exp(XN(i)), 0.0);
		}
		for (int i = num_path - 1; i >= 0; --i) {
			XN(i) = log(S0) + (num_path - i) * dx;
			FN(i) = max(-K + exp(XN(i)), 0.0);
		}
	}
	MatrixXd A = MatrixXd::Zero(total_num_path, total_num_path);
	
	A(0,0) = pu;
	A(0,1) = pm;
	A(0,2) = pd;
	A(total_num_path-1,total_num_path-3) = pu;
	A(total_num_path-1,total_num_path-2) = pm;
	A(total_num_path-1,total_num_path-1) = pd;
	for (int i = 1; i < total_num_path-1; ++i) {
		A(i,i - 1) = pu;
		A(i,i) = pm;
		A(i,i + 1) = pd;
	}
	VectorXd Fi;
	VectorXd temp;
	Fi = FN;
	VectorXd B = VectorXd::Zero(total_num_path);
	if (call == 0) {
		B(total_num_path - 1) = -exp(XN(total_num_path - 1)) + exp(XN(total_num_path - 2));
	}
	else {
		B(0) = exp(XN(0)) - exp(XN(1));
	}
	if (euro == 1) {
		// european option
		for (int i = 0; i < N; ++i) {
			temp = A * Fi + B;
			Fi = temp;
		}
	}
	else {
		// american option
		for (int i = 0; i < N; ++i) {
			temp = A * Fi + B;
			if (call == 1) {
				for (int j = 0; j < total_num_path; ++j) {
					temp(j) = max(temp(j), exp(XN(j))-K);
				}
			}
			else {
				for (int j = 0; j < total_num_path; ++j) {
					temp(j) = max(temp(j), -exp(XN(j)) + K);
				}
			}
			Fi = temp;
		}
	}
	return Fi(num_path);
}

double IFDM(double S0, int k, int euro, int call) {
	double sigma = 0.2;
	double T = 0.5;
	double K = 10;
	double dt = 0.002;
	double N = T / dt;
	double dx = sigma * sqrt(k*dt);
	double r = 0.04;
	double Smax = 16;
	double Smin = 4;
	double pu = -0.5* dt * (sigma*sigma / (dx*dx) + (r - sigma * sigma / 2) / (dx));
	double pm = 1 + dt * sigma*sigma / (dx*dx) + r * dt;
	double pd = -0.5*dt * (sigma*sigma / (dx*dx) - (r - sigma * sigma / 2) / (dx));
	int num_path = (log(S0) - log(Smin)) / dx;
	int total_num_path = 2 * num_path + 1;
	VectorXd XN = VectorXd::Zero(total_num_path);
	VectorXd BN(total_num_path);

	if (call == 0) {
		for (int i = num_path; i < 2 * num_path + 1; ++i) {
			XN(i) = log(S0) - (i - num_path) * dx;
			BN(i) = max(K - exp(XN(i)), 0.0);
		}
		for (int i = num_path - 1; i >= 0; --i) {
			XN(i) = log(S0) + (num_path - i) * dx;
			BN(i) = max(K - exp(XN(i)), 0.0);
		}
		BN(0) = 0;
		BN(total_num_path- 1) = -exp(XN(total_num_path - 1)) + exp(XN(total_num_path - 2));
	}
	else {
		for (int i = num_path; i < 2 * num_path + 1; ++i) {
			XN(i) = log(S0) - (i - num_path) * dx;
			BN(i) = max(-K + exp(XN(i)), 0.0);
		}
		for (int i = num_path - 1; i >= 0; --i) {
			XN(i) = log(S0) + (num_path - i) * dx;
			BN(i) = max(-K + exp(XN(i)), 0.0);
		}
		BN(0) = exp(XN(0)) - exp(XN(1));
		BN(total_num_path - 1) = 0;
	}

	MatrixXd A = MatrixXd::Zero(total_num_path, total_num_path);
	A(0,0) = 1;
	A(0,1) = -1;
	A(total_num_path - 1,total_num_path - 2) = 1;
	A(total_num_path - 1,total_num_path - 1) = -1;
	for (int i = 1; i < total_num_path - 1; ++i) {
		A(i,i - 1) = pu;
		A(i,i) = pm;
		A(i,i + 1) = pd;
	}
	
	VectorXd temp;
	MatrixXd Ainv = A.inverse();
	VectorXd Fi = Ainv * BN;
	if (euro == 1) {
		for (int i = 0; i < N; ++i) {
			temp = Fi;
			if (call == 0) {
				temp(0) = 0;
				temp(total_num_path - 1) = -exp(XN(total_num_path - 1)) + exp(XN(total_num_path - 2));
			}
			else {
				temp(0) = exp(XN(0)) - exp(XN(1));
				temp(total_num_path - 1) = 0;
			}
			Fi = Ainv * temp;
		}
	}
	else {
		for (int i = 0; i < N; ++i) {
			temp = Fi;
			
			if (call == 0) {
				for (int j = 0; j < total_num_path; ++j) {
					temp(j) = max(temp(j), -exp(XN(j)) + K);
				}
				temp(0) = 0;
				temp(total_num_path - 1) = -exp(XN(total_num_path - 1)) + exp(XN(total_num_path - 2));
			}
			else {
				for (int j = 0; j < total_num_path; ++j) {
					temp(j) = max(temp(j), exp(XN(j)) - K);
				}
				temp(0) = exp(XN(0)) - exp(XN(1));
				temp(total_num_path - 1) = 0;
			}
			Fi = Ainv * temp;
		}
	}
	return Fi(num_path);
}

double CNFDM(double S0, int k) {
	double sigma = 0.2;
	double T = 0.5;
	double K = 10;
	double dt = 0.002;
	double N = T / dt;
	double dx = sigma * sqrt(k*dt);
	double r = 0.04;
	double Smax = 16;
	double Smin = 4;
	double pu = -0.25* dt * (sigma*sigma / (dx*dx) + (r - sigma * sigma / 2) / (dx));
	double pm = 1 + dt * sigma*sigma / (2*dx*dx) + r * dt/2;
	double pd = -0.25*dt * (sigma*sigma / (dx*dx) - (r - sigma * sigma / 2) / (dx));
	int num_path = (log(S0) - log(Smin)) / dx;
	int total_num_path = 2 * num_path + 1;
	VectorXd XN = VectorXd::Zero(total_num_path);
	VectorXd CN(total_num_path);
	for (int i = num_path; i < 2 * num_path + 1; ++i) {
		XN(i) = log(S0) - (i - num_path) * dx;
		CN(i) = max(K - exp(XN(i)), 0.0);
	}
	for (int i = num_path - 1; i >= 0; --i) {
		XN(i) = log(S0) + (num_path - i) * dx;
		CN(i) = max(K - exp(XN(i)), 0.0);
	}

	MatrixXd Az = MatrixXd::Zero(total_num_path, total_num_path);
	Az(0, 0) = 1;
	Az(0, 1) = -1;
	Az(total_num_path - 1, total_num_path - 2) = 1;
	Az(total_num_path - 1, total_num_path - 1) = -1;
	MatrixXd A = MatrixXd::Zero(total_num_path, total_num_path);
	A(0, 0) = 1;
	A(0, 1) = -1;
	A(total_num_path - 1, total_num_path - 2) = 1;
	A(total_num_path - 1, total_num_path - 1) = -1;
	for (int i = 1; i < total_num_path - 1; ++i) {
		A(i, i - 1) = pu;
		A(i, i) = pm;
		A(i, i + 1) = pd;

		Az(i, i - 1) = -pu;
		Az(i, i) = -(pm-2);
		Az(i, i + 1) = -pd;
	}
	VectorXd Zi = Az * CN;
	
	VectorXd temp;
	MatrixXd Ainv = A.inverse();
	VectorXd Ci;
	for (int i = 0; i < N; ++i) {
		temp = Zi;
		temp(0) = 0;
		temp(total_num_path - 1) = -exp(XN(total_num_path - 1)) + exp(XN(total_num_path - 2));
		Ci = Ainv * temp;
		Zi = Az * Ci;
	}
	return Ci(num_path);
}

double EFDM_S(double S0, double dS, double alpha, int call) {
	double sigma = 0.2;
	double T = 0.5;
	double K = 10;
	double dt = 0.002;
	double N = T / dt;
	
	double r = 0.04;
	double Smax = 16;
	double Smin = 4;

	int num_path = S0/ dS;
	int total_num_path = 2 * num_path + 1;
	VectorXd SN = VectorXd::Zero(total_num_path);
	VectorXd FN(total_num_path);
	if (call == 0) {
		// put option
		for (int i = 0; i < total_num_path; ++i) {
			SN(i) = (total_num_path - 1 -i) * dS ;
			FN(i) = max(K - SN(i), 0.0);
		}
	}
	else {
		// call option
		for (int i = 0; i < total_num_path; ++i) {
			SN(i) = (total_num_path - 1 - i) * dS;
			FN(i) = max(-K + SN(i), 0.0);
		}
	}
	
	MatrixXd A = MatrixXd::Zero(total_num_path, total_num_path);
	MatrixXd Ad = MatrixXd::Zero(total_num_path, total_num_path);
	Ad(0, 0) = 1;
	Ad(0, 1) = -1;
	Ad(total_num_path - 1, total_num_path - 2) = 1;
	Ad(total_num_path - 1, total_num_path - 1) = -1;
	A(0, 0) = 1;
	A(0, 1) = -1;
	A(total_num_path - 1, total_num_path - 2) = 1;
	A(total_num_path - 1, total_num_path - 1) = -1;
	for (int i = 1; i < total_num_path - 1; ++i) {
		double j = total_num_path - 1 - i;
		//a3
		A(i, i - 1) = r * j*(1 - alpha) / 2 + sigma * sigma*j*j*(1 - alpha) / 2;
		//a2
		A(i, i) = -1 / dt - sigma * sigma*j*j*(1 - alpha) - r * (1 - alpha);
		//a1
		A(i, i + 1) = -r * j*(1 - alpha) / 2 + sigma * sigma*j*j*(1 - alpha) / 2;
		//b3
		Ad(i, i - 1) = -(r * j*(alpha) / 2 + sigma * sigma*j*j*(alpha) / 2);
		//b2
		Ad(i, i) = -(1 / dt - sigma * sigma*j*j*(alpha) - r * (alpha));
		//b1
		Ad(i, i + 1) = -(-r * j*(alpha) / 2 + sigma * sigma*j*j*(alpha) / 2);

	}

	VectorXd Di = Ad * FN;
	
	VectorXd temp;
	MatrixXd Ainv = A.inverse();
	VectorXd Ci;
	for (int i = 0; i < N; ++i) {
		temp = Di;
		if (call == 0) {
			temp(0) = 0;
			temp(total_num_path - 1) = -SN(total_num_path - 1) + SN(total_num_path - 2);
			Ci = Ainv * temp;
			for (int j = 0; j < total_num_path; ++j) {
				Ci(j) = max(Ci(j), K - SN(j));
			}
		}
		else {
			temp(0) = SN(0) - SN(1);
			temp(total_num_path - 1) = 0;
			Ci = Ainv * temp;
			for (int j = 0; j < total_num_path; ++j) {
				Ci(j) = max(Ci(j), -K + SN(j));
			}
		}
		//cout << Ci << endl;
		Di = Ad * Ci;
	}
	return Ci(num_path);

}