#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "mbs.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>

using namespace std;

const double pi = atan(1) * 4;

double MBS(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double x, int64_t seed) {
	double r = WAC / 12;
	int N = T * 12;
	double dt = 1.0/12;
	int sim_num = 10000;
	double month[] = { 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98 };
	double* out = new double[sim_num];
	for (int i = 0; i < sim_num; ++i) {
		double* rt = new double[N + 121];
		double* uni = getLGM(N+120, seed + 1000 * i);
		double* rand = box_muller(uni, N+120);
		double* d_t = new double[N];
		double* c_t = new double[N];
		double* PVt_1 = new double[N + 1];
		double* TPP_t = new double[N];
		double* CPR_t = new double[N];
		double* IP_t = new double[N];
		double* rt_1_10 = new double[N];
		double* RI_t = new double[N];
		double* BU_t = new double[N];
		double* SG_t = new double[N];
		double* SY_t = new double[N];
		PVt_1[0] = PV0;
		rt[0] = r0;
		for (int j = 1; j < N + 121;++j) {
			rt[j] = rt[j - 1] + kappa * (r_bar - rt[j - 1])*dt + sigma * sqrt(abs(rt[j - 1])) * sqrt(dt) * rand[j - 1];
		}
		
		double sum = 0;
		for (int j = 0; j < N; ++j) {
			//t starts from 1
			d_t[j] = exp(-trapzoid_method(0, j + 1, rt, dt));
			rt_1_10[j] = -1.0 / 10 * trapzoid_method(j, j+120, rt, dt);
			RI_t[j] = 0.28 + 0.14 * atan(-8.57 + 430 * (WAC - rt_1_10[j]));
			BU_t[j] = 0.3 + 0.7* PVt_1[j] / PV0;
			SG_t[j] = min(1.0, (double)(j + 1) / 30);
			SY_t[j] = month[j % 12];
			CPR_t[j] = RI_t[j] * BU_t[j] * SG_t[j] * SY_t[j];
			TPP_t[j] = PVt_1[j] * r * (1.0 / (1 - pow((1 + r), (-N + j))) - 1) + (PVt_1[j] - PVt_1[j] * r*(1.0 / (1 - pow((1 + r), (-N + j))) - 1))*(1.0 - pow(1 - CPR_t[j], 1.0 / 12));
			IP_t[j] = PVt_1[j] * r;
			PVt_1[j + 1] = PVt_1[j] - TPP_t[j];
			c_t[j] = TPP_t[j] + IP_t[j];
			sum = sum + d_t[j] * c_t[j] * exp(-x*(j+1)*dt);
			
		}
		out[i] = sum;
		delete[] rt ;
		delete[] uni 	;
		delete[] rand 	;
		delete[] d_t 	;
		delete[] c_t 	;
		delete[] PVt_1	;
		delete[] TPP_t	;
		delete[] CPR_t	;
		delete[] IP_t 	;
		delete[] rt_1_10	;
		delete[] RI_t 	;
		delete[] BU_t 	;
		delete[] SG_t 	;
		delete[] SY_t 	;
	}
	double result = calc_mean(out, sim_num);
	delete[] out;
	return result;
}

double OAS(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double mkt_p, int64_t seed) {
	
	double r = WAC / 12;
	int N = T * 12;
	double dt = 1.0 / 12;
	int sim_num = 10000;
	double month[] = { 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98 };
	double** all_ct = new double*[sim_num];
	double** all_dt = new double*[sim_num];
	//double** out = new double*[sim_num];
	double* xt = new double[N];
	for (int i = 0; i < sim_num; ++i) {
		//out[i] = new double[N];
		all_ct[i] = new double[N];
		all_dt[i] = new double[N];
		double* rt = new double[N + 121];
		double* uni = getLGM(N + 120, seed + 1000 * i);
		double* rand = box_muller(uni, N + 120);
		double* d_t = new double[N];
		double* c_t = new double[N];
		double* PVt_1 = new double[N + 1];
		double* TPP_t = new double[N];
		double* CPR_t = new double[N];
		double* IP_t = new double[N];
		double* rt_1_10 = new double[N];
		double* RI_t = new double[N];
		double* BU_t = new double[N];
		double* SG_t = new double[N];
		double* SY_t = new double[N];
		PVt_1[0] = PV0;
		rt[0] = r0;
		for (int j = 1; j < N + 121;++j) {
			rt[j] = rt[j - 1] + kappa * (r_bar - rt[j - 1])*dt + sigma * sqrt(abs(rt[j - 1])) * sqrt(dt) * rand[j - 1];
		}

		double sum = 0;
		for (int j = 0; j < N; ++j) {
			//t starts from 1
			all_dt[i][j] = exp(-trapzoid_method(0, j + 1, rt, dt));
			rt_1_10[j] = -1.0 / 10 * trapzoid_method(j, j + 120, rt, dt);
			RI_t[j] = 0.28 + 0.14 * atan(-8.57 + 430 * (WAC - rt_1_10[j]));
			BU_t[j] = 0.3 + 0.7* PVt_1[j] / PV0;
			SG_t[j] = min(1.0, (double)(j + 1) / 30);
			SY_t[j] = month[j % 12];
			CPR_t[j] = RI_t[j] * BU_t[j] * SG_t[j] * SY_t[j];
			TPP_t[j] = PVt_1[j] * r * (1.0 / (1 - pow((1 + r), (-N + j))) - 1) + (PVt_1[j] - PVt_1[j] * r*(1.0 / (1 - pow((1 + r), (-N + j))) - 1))*(1 - pow(1 - CPR_t[j], 1.0 / 12));
			IP_t[j] = PVt_1[j] * r;
			PVt_1[j + 1] = PVt_1[j] - TPP_t[j];
			all_ct[i][j] = TPP_t[j] + IP_t[j];
			//out[i][j] = all_dt[i][j] * all_ct[i][j];

			//sum = sum + d_t[j] * c_t[j];

		}
		//out[i] = sum;
		delete[] rt;
		delete[] uni;
		delete[] rand;
		delete[] d_t;
		delete[] c_t;
		delete[] PVt_1;
		delete[] TPP_t;
		delete[] CPR_t;
		delete[] IP_t;
		delete[] rt_1_10;
		delete[] RI_t;
		delete[] BU_t;
		delete[] SG_t;
		delete[] SY_t;
	}

	double result = 0.0;
	double left = -1.0;
	double right = 0.0;
	double x = 0.0;
	//double x = (left + right) / 2;
	double* out = new double[sim_num];
	int count = 0;
	while (abs(result - mkt_p) > 0.00000001) {
		++count;

		if (result < mkt_p) {
			right = x;
			x = (left + right) / 2;
		}
		else {
			left = x;
			x = (left + right) / 2;
		}
		for (int i = 0; i < N; ++i) {
			xt[i] = exp(-(i + 1)*x *dt);
		}
		for (int i = 0; i < sim_num; ++i) {
			double sum = 0;
			for (int j = 0; j < N; ++j) {
				sum = sum + all_dt[i][j] * all_ct[i][j] * xt[j];
			}
			out[i] = sum;
		}
		result = calc_mean(out, sim_num);
		
	}
	return x;
}


double MBS_explicit(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double x) {
	double r = WAC / 12;
	int N = T * 12;
	double dt = 1.0 / 12;
	double month[] = { 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98 };
	double* rt = new double[N + 121];
	double* B_0t = new double[N];
	double* A_0t = new double[N];
	double* d_t = new double[N];
	double* c_t = new double[N];
	double* PVt_1 = new double[N + 1];
	double* TPP_t = new double[N];
	double* CPR_t = new double[N];
	double* IP_t = new double[N];
	double* rt_1_10 = new double[N];
	double* RI_t = new double[N];
	double* BU_t = new double[N];
	double* SG_t = new double[N];
	double* SY_t = new double[N];
	PVt_1[0] = PV0;
	rt[0] = r0;
	for (int i = 1; i < N + 121; ++i) {
		rt[i] = rt[0] * exp(-kappa * (i)*dt) + r_bar * (1 - exp(-kappa * i*dt));
		//cout << rt[i] << endl;
	}
	//cout << trapzoid_method(0,N,rt, dt) << endl;
	double sum = 0;
	double h1 = sqrt(kappa*kappa + 2 * sigma*sigma);
	double h2 = (kappa + h1) / 2;
	double h3 = 2 * kappa*r_bar / (sigma*sigma);
	for (int j = 0; j < N; ++j) {
		//t starts from 1
		B_0t[j] = (exp(h1*(j + 1)*dt) - 1) / (h2*(exp(h1*(j + 1)*dt) - 1) + h1);
		A_0t[j] = pow((h1* exp(h2*(j + 1)*dt)) / (h2*(exp(h1*(j + 1)*dt) - 1) + h1),h3);
		d_t[j] = A_0t[j] * exp(-B_0t[j]*r0);
		double B_j10 = (exp(h1 * 10) - 1) / (h2*(exp(h1 * 10) - 1) + h1);
		double A_j10 = pow((h1*exp(h2 * 10)) / (h2*(exp(h1 * 10) - 1) + h1), h3);
		rt_1_10[j] = 1.0 / 10 * log(A_j10 * exp(B_j10*rt[j]));
		RI_t[j] = 0.28 + 0.14 * atan(-8.57 + 430 * (WAC - rt_1_10[j]));
		BU_t[j] = 0.3 + 0.7* PVt_1[j] / PV0;
		SG_t[j] = min(1.0, (double)(j + 1) / 30);
		SY_t[j] = month[j % 12];
		CPR_t[j] = RI_t[j] * BU_t[j] * SG_t[j] * SY_t[j];
		TPP_t[j] = PVt_1[j] * r * (1.0 / (1 - pow((1 + r), (-N + j))) - 1) + (PVt_1[j] - PVt_1[j] * r*(1.0 / (1 - pow((1 + r), (-N + j))) - 1))*(1.0 - pow(1 - CPR_t[j], 1.0 / 12));
		IP_t[j] = PVt_1[j] * r;
		PVt_1[j + 1] = PVt_1[j] - TPP_t[j];
		c_t[j] = TPP_t[j] + IP_t[j];
		sum = sum + d_t[j] * c_t[j] * exp(-x * (j + 1)*dt);
	}
	
	delete[] rt;
	delete[] d_t;
	delete[] c_t;
	delete[] PVt_1;
	delete[] TPP_t;
	delete[] CPR_t;
	delete[] IP_t;
	delete[] rt_1_10;
	delete[] RI_t;
	delete[] BU_t;
	delete[] SG_t;
	delete[] SY_t;
	delete[] B_0t;
	delete[] A_0t;
	return sum;
}

double OAS_explicit(double WAC, int T, double PV0, double r0, double kappa, double r_bar, double sigma, double mkt_p) {
	double r = WAC / 12;
	int N = T * 12;
	double dt = 1.0 / 12;
	double month[] = { 0.94, 0.76, 0.74, 0.95, 0.98, 0.92, 0.98, 1.10, 1.18, 1.22, 1.23, 0.98 };
	double* rt = new double[N + 121];
	double* B_0t = new double[N];
	double* A_0t = new double[N];
	double* d_t = new double[N];
	double* c_t = new double[N];
	double* PVt_1 = new double[N + 1];
	double* TPP_t = new double[N];
	double* CPR_t = new double[N];
	double* IP_t = new double[N];
	double* rt_1_10 = new double[N];
	double* RI_t = new double[N];
	double* BU_t = new double[N];
	double* SG_t = new double[N];
	double* SY_t = new double[N];
	PVt_1[0] = PV0;
	rt[0] = r0;
	for (int i = 1; i < N + 121; ++i) {
		rt[i] = rt[0] * exp(-kappa * (i)*dt) + r_bar * (1 - exp(-kappa * i*dt));
	}

	double sum = 0;
	double h1 = sqrt(kappa*kappa + 2 * sigma*sigma);
	double h2 = (kappa + h1) / 2;
	double h3 = 2 * kappa*r_bar / (sigma*sigma);
	for (int j = 0; j < N; ++j) {
		//t starts from 1
		B_0t[j] = (exp(h1*(j + 1)*dt) - 1) / (h2*(exp(h1*(j + 1)*dt) - 1) + h1);
		A_0t[j] = pow((h1* exp(h2*(j + 1)*dt)) / (h2*(exp(h1*(j + 1)*dt) - 1) + h1), h3);
		d_t[j] = A_0t[j] * exp(-B_0t[j] * r0);
		double B_j10 = (exp(h1 * 10) - 1) / (h2*(exp(h1 * 10) - 1) + h1);
		double A_j10 = pow((h1*exp(h2 * 10)) / (h2*(exp(h1 * 10) - 1) + h1), h3);
		rt_1_10[j] = 1.0 / 10 * log(A_j10 * exp(B_j10*rt[j]));
		RI_t[j] = 0.28 + 0.14 * atan(-8.57 + 430 * (WAC - rt_1_10[j]));
		BU_t[j] = 0.3 + 0.7* PVt_1[j] / PV0;
		SG_t[j] = min(1.0, (double)(j + 1) / 30);
		SY_t[j] = month[j % 12];
		CPR_t[j] = RI_t[j] * BU_t[j] * SG_t[j] * SY_t[j];
		TPP_t[j] = PVt_1[j] * r * (1.0 / (1 - pow((1 + r), (-N + j))) - 1) + (PVt_1[j] - PVt_1[j] * r*(1.0 / (1 - pow((1 + r), (-N + j))) - 1))*(1.0 - pow(1 - CPR_t[j], 1.0 / 12));
		IP_t[j] = PVt_1[j] * r;
		PVt_1[j + 1] = PVt_1[j] - TPP_t[j];
		c_t[j] = TPP_t[j] + IP_t[j];
		sum = sum + d_t[j] * c_t[j];

	}

	delete[] rt;
	delete[] PVt_1;
	delete[] TPP_t;
	delete[] CPR_t;
	delete[] IP_t;
	delete[] rt_1_10;
	delete[] RI_t;
	delete[] BU_t;
	delete[] SG_t;
	delete[] SY_t;
	delete[] B_0t;
	delete[] A_0t;

	double* xt = new double[N];
	
	double left = -1.0;
	double right = 0.0;
	double x = 0.0;
	//double x = (left + right) / 2;
	int count = 0;
	while (abs(sum - mkt_p) > 0.00000001) {
		++count;
		
		if (sum < mkt_p) {
			right = x;
			x = (left + right) / 2;
		}
		else {
			left = x;
			x = (left + right) / 2;
		}
		sum = 0;
		for (int i = 0; i < N; ++i) {
			xt[i] = exp(-(i + 1)*x *dt);
			sum = sum + d_t[i] * c_t[i] * xt[i];
		}
		

	}
	delete[] xt;
	delete[] c_t;
	delete[] d_t;
	return x;
}

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

double trapzoid_method(int a, int b, double* integrand, double dt) {
	double sum = 0;
	int i = a;
	while (i < b) {
		sum = sum + (integrand[i] + integrand[i + 1])*dt / 2;
		++i;
	}
	return sum;
}

double calc_mean(double *nums, int n) {
	double sum_of_mean = 0;
	for (int i = 0;i < n;++i) {
		sum_of_mean = sum_of_mean + nums[i];
	}
	double mean = sum_of_mean / n;
	return mean;
}