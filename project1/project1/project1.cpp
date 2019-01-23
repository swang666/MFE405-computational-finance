// project1.cpp 
// author: Sumeng Wang

#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include <ctime>
#include <string>
#include <array>
#include <chrono> 
#include <fstream>
using namespace std::chrono;
using namespace std;
double * getLGM(int n);
double sum(double *nums, int n);
double * gen_bern(double p, int n);
double * gen_binom(double p, int n, int c);
double * box_muller(double *nums, int n);
double * polar_mar(double *nums, int n, int* length);
double * polar_mar2(double *nums, int n);

double fact(int n)
{
	double res = 1;
	for (int i = 2; i <= n; i++)
		res = res * i;
	return res;
}

double * getLGM(int n) {
	int64_t m = pow(2, 31) - 1;
	int a = pow(7, 5);
	int b = 0;
	int64_t X0 = 12345; //some arbitrary number
	double* r = new double[n];
	int64_t* x = new int64_t[n];
	x[0] = X0;
	r[0] = (double)(X0+0.5)/m;
	for (int i = 1; i < n; ++i) {
		x[i] = (a * x[i-1]) % m;
		r[i] = (double)(x[i]+0.5) / m;
	}
	return r;
}

double sum(double *nums, int n) {
	double sum_of_mean = 0;
	for (int i = 0;i < n;++i) {
		sum_of_mean = sum_of_mean + nums[i];
	}
	return sum_of_mean;
}

double calc_mean(double *nums, int n) {
	double sum_of_mean = 0;
	for (int i = 0;i < n;++i) {
		sum_of_mean = sum_of_mean + nums[i];
	}
	double mean = sum_of_mean / n;
	return mean;
}

double calc_std(double *nums, int n) {
	double sum_of_std = 0;
	double mean = calc_mean(nums, n);
	for (int i = 0;i < n;++i) {
		sum_of_std = sum_of_std + pow((nums[i] - mean), 2.0);
	}
	double sigma = sqrt(sum_of_std / n);
	return sigma;
}

double * gen_bern(double p, int n) {
	double *rand_num = getLGM(n);
	double* x = new double[n];
	for (int i = 0; i < n; ++i) {
		if (rand_num[i] <= p) {
			x[i] = 1;
		}
		else {
			x[i] = 0;
		}
	}
	return x;
}

double * gen_binom(double p, int n, int c) {
	double* binom = new double[c];
	double *rand_bern = gen_bern(p, n*c);
	double sum_bern;
	for (int i = 0; i < c; ++i) {
		sum_bern = 0;
		for (int j = 0; j < n;++j) {
			sum_bern = sum_bern + rand_bern[i*n + j];
		}
		binom[i] = sum_bern;
	}
	free(rand_bern);
	return binom;
}

double * gen_exp(double lambda, int n) {
	double *rand_num = getLGM(n);
	double *rand_exp = new double[n];
	for (int i = 0; i < n; ++i) {
		rand_exp[i] = -lambda * log(rand_num[i]);
	}
	return rand_exp;
}

int main()
{
	//Question 1
	//a
	cout << "Question 1" << endl;
	double *randnum = getLGM(10000);
	double sum_of_std = 0;
	double mean = calc_mean(randnum, 10000);
	double sigma = calc_std(randnum, 10000);
	cout << "a" << endl;
	cout << "mean: " << mean << endl;
	cout << "stdev: " << sigma << endl;

	//b
	srand(12345);
	double randnum_b[10000];
	for (int i = 0;i < 10000;++i) {
		randnum_b[i] = rand() / (RAND_MAX + 1.);
	}
	double mean_b = calc_mean(randnum_b, 10000);
	double sigma_b = calc_std(randnum_b, 10000);
	cout << "b" << endl;
	cout << "mean: " << mean_b << endl;
	cout << "stdev: " << sigma_b << endl << endl;

	//Question 2
	cout << "Question 2" << endl;
	//a
	double x2[10000];
	ofstream myfile;
	myfile.open("ques2.csv", fstream::out);
	int occurence[4] = { 0, 0, 0, 0 };
	for (int i = 0;i < 10000; ++i) {
		if (randnum[i] < 0.3) {
			x2[i] = -1.0;
			++occurence[0];
		}
		else if (randnum[i] < 0.65) {
			x2[i] = 0.0;
			++occurence[1];
		}
		else if (randnum[i] < 0.85) {
			x2[i] = 1.0;
			++occurence[2];
		}
		else {
			x2[i] = 2.0;
			++occurence[3];
		}
		myfile << to_string(x2[i]) << '\n';
	}
	myfile.close();
	//b
	cout << "b" << endl;
	double mean_x2 = calc_mean(x2, 10000);
	double sigma_x2 = calc_std(x2, 10000);
	cout << "mean: " << mean_x2 << endl;
	cout << "stdev: " << sigma_x2 << endl << endl;

	//Question 3
	cout << "Question 3" << endl;
	double* binom = gen_binom(0.64, 44, 1000);
	myfile.open("ques3.csv", fstream::out);
	int num_success = 0;
	for (int i = 0; i < 1000; ++i) {
		if (binom[i] >= 40) {
			++num_success;
		}
		myfile << to_string(binom[i]) << '\n';
	}
	myfile.close();
	cout << "probability of the sample is: " << (double) num_success/1000 << endl;
	double prob = fact(44) / fact(4) / fact(40) * pow(0.64, 40) * pow(0.36, 4) +
		fact(44) / fact(3) / fact(41) * pow(0.64, 41) * pow(0.36, 3) +
		fact(44) / fact(2) / fact(42) * pow(0.64, 42) * pow(0.36, 2) +
		fact(44) / fact(43) * pow(0.64, 43) * pow(0.36, 1) +
		fact(44) / fact(44) * pow(0.64, 44) * pow(0.36, 0);
	cout << "theoretical probability is: " << prob << endl << endl;

	//Question 4
	cout << "Question 4" << endl;
	//a
	double* exp = gen_exp(1.5, 10000);
	myfile.open("ques4.csv", fstream::out);
	cout << "b" << endl;
	//b
	int num_greathan_1 = 0;
	for (int i = 0; i < 10000; ++i) {
		if (exp[i] >= 1) {
			++num_greathan_1;
		}
		myfile << to_string(exp[i]) << '\n';
	}
	myfile.close();

	double prob_1 = (double)num_greathan_1 / 10000;
	cout << "P(X >= 1) = " << prob_1 << endl;

	int num_greathan_4 = 0;
	for (int i = 0; i < 10000; ++i) {
		if (exp[i] >= 4) {
			++num_greathan_4;
		}
	}

	double prob_4 = (double)num_greathan_4 / 10000;
	cout << "P(X >= 4) = " << prob_4 << endl;

	//c
	cout << "c" << endl;
	double exp_mean = calc_mean(exp, 10000);
	double exp_sd = calc_std(exp, 10000);
	cout << "mean: " << exp_mean << endl;
	cout << "stdev: " << exp_sd << endl << endl;

	//Question 5
	cout << "Question 5" << endl;
	//a
	double *rand_num_5 = getLGM(5000);

	//b
	double *std_norm_box = box_muller(rand_num_5, 5000);

	//c
	cout << "c" << endl;
	double box_mean = calc_mean(std_norm_box, 5000);
	double box_sd = calc_std(std_norm_box, 5000);
	cout << "box-muller mean: " << box_mean << endl;
	cout << "box-muller stdev: " << box_sd << endl;

	//d
	int length;
	double *temp = polar_mar(rand_num_5, 5000, &length);
	double *std_norm_polar = new double[length];
	for (int i = 0; i < length; ++i) {
		std_norm_polar[i] = temp[i];
	}

	//e
	cout << "e" << endl;
	double polar_mean = calc_mean(std_norm_polar, length);
	double polar_sd = calc_std(std_norm_polar, length);
	cout << "polar-marsaglia mean: " << polar_mean << endl;
	cout << "polar-marsaglia stdev: " << polar_sd << endl;

	//f
	cout << "f" << endl;
	int k = 0;
	for (int i = 0; i < 1000; ++i) {
		double *rand_num = getLGM(20000);
		auto start = high_resolution_clock::now();
		double *std_norm_box = box_muller(rand_num, 5000);
		// Get ending timepoint 
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		auto start2 = high_resolution_clock::now();
		double *std_norm_polar = polar_mar2(rand_num, 5000);
		auto stop2 = high_resolution_clock::now();
		auto duration2 = duration_cast<microseconds>(stop2 - start2);
		if (duration.count() > duration2.count()) {
			// when polar is faster
			++k;
		}
	}

	cout << (double)k / 1000 * 100 << "% of times, polar-marsaglia is faster" << endl;

}

double * polar_mar2(double *nums, int n) {
	double z1, z2, v1, v2, w;
	double *temp = new double[n];
	int count = 0;
	int i = 0;
	while (count < n) {
		v1 = 2 * nums[i] - 1;
		v2 = 2 * nums[i + 1] - 1;
		w = v1 * v1 + v2 * v2;
		if (w <= 1) {
			z1 = v1 * sqrt(-2 * log(w) / w);
			z2 = v2 * sqrt(-2 * log(w) / w);
			temp[count] = z1;
			temp[count + 1] = z2;
			count = count + 2;
		}
		i = i + 2;
	}
	return temp;
}

double * polar_mar(double *nums, int n, int* length) {
	double z1, z2, v1, v2, w;
	double *temp = new double[n];
	int count = 0;
	for (int i = 0; i < n; i = i + 2) {
		v1 = 2 * nums[i] - 1;
		v2 = 2 * nums[i + 1] - 1;
		w = pow(v1, 2) + pow(v2, 2);
		if (w <= 1) {
			z1 = v1 * sqrt(-2 * log(w) / w);
			z2 = v2 * sqrt(-2 * log(w) / w);
			temp[count] = z1;
			temp[count + 1] = z2;
			count = count + 2;
		}
	}
	*length = count - 2;
	return temp;
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