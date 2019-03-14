// Project9.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "mbs.h"
#include <fstream>
#include <string>

using namespace std;

int main()
{
	//ofstream myfile;
	cout << MBS_explicit(0.08, 30, 100000, 0.078, 0.6, 0.08, 0.12, 0) << endl;
	
	double kappas[] = { 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
	double* q1b = new double[7];
	double rbars[] = { 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09 };
	double* q1c = new double[7];
	double sigmas[] = { 0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2 };
	double* q1d = new double[11];
	//myfile.open("q1bc_ex.csv", fstream::out);
	//myfile << "kappa" << ',' << "price1" << ',' << "r_bar" << ',' << "price2" << endl;
	for (int i = 0; i < 7; ++i) {
		q1b[i] = MBS_explicit(0.08, 30, 100000, 0.078, kappas[i], 0.08, 0.12, 0);
		q1c[i] = MBS_explicit(0.08, 30, 100000, 0.078, 0.6, rbars[i], 0.12, 0);
		//myfile << to_string(kappas[i]) << ',' << to_string(q1b[i]) <<','<< to_string(rbars[i]) << ','<< to_string(q1c[i]) << endl;
	}
	//myfile.close();
	//myfile.open("q1d_ex.csv", fstream::out);
	//myfile << "sigma" << ',' << "price" << endl;
	for (int i = 0; i < 11;++i) {
		q1d[i] = MBS_explicit(0.08, 30, 100000, 0.078, 0.6, 0.08, sigmas[i], 0);
		//myfile << to_string(sigmas[i]) << ',' << to_string(q1d[i]) << endl;
	}
	//myfile.close();

	double x = OAS(0.08, 30, 100000, 0.078, 0.6, 0.08, 0.12, 110000,12345);
	cout << "OAS: " << x << endl;

	double y = 0.0005;
	double y_plus = x + y;
	double y_minus = x - y;
	double P_plus = MBS(0.08, 30, 100000, 0.078, 0.6, 0.08, 0.12, y_plus,12345);
	double P_minus = MBS(0.08, 30, 100000, 0.078, 0.6, 0.08, 0.12, y_minus,12345);
	cout << "P+: " << P_plus << endl;
	cout << "P-: " << P_minus << endl;
	double D = (P_minus - P_plus) / (2 * y * 110000);
	double C = (P_plus + P_minus - 2 * 110000) / (2 * 110000 * y * y);
	cout << "Duration: " << D << endl;
	cout << "Convexity: " << C << endl;


	double x2 = OAS_explicit(0.08, 30, 100000, 0.078, 0.6, 0.08, 0.12, 110000);
	cout << "OAS: " << x2 << endl;

	double y2 = 0.0005;
	double y_plus2 = x2 + y2;
	double y_minus2 = x2 - y2;
	double P_plus2 = MBS_explicit(0.08, 30, 100000, 0.078, 0.6, 0.08, 0.12, y_plus2);
	double P_minus2 = MBS_explicit(0.08, 30, 100000, 0.078, 0.6, 0.08, 0.12, y_minus2);
	cout << "P+: " << P_plus2 << endl;
	cout << "P-: " << P_minus2 << endl;
	double D2 = (P_minus2 - P_plus2) / (2 * y2 * 110000);
	double C2 = (P_plus2 + P_minus2 - 2 * 110000) / (2 * 110000 * y2 * y2);
	cout << "Duration: " << D2 << endl;
	cout << "Convexity: " << C2 << endl;
}
