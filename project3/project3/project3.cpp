// project3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "rng.h"
#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;

int main()
{
	int64_t seed = 12345;
	cout << "Question 1" << endl;

	double * ques1 = question1(seed);
	cout << "Prob: " << ques1[0] << endl;
	cout << "E1: " << ques1[1] << endl;
	cout << "E2: " << ques1[2] << endl;
	cout << "E3: " << ques1[3] << endl << endl;
	cout << "Question 2" << endl;
	double* ques2 = question2(seed);
	cout << "E1: " << ques2[0] << endl;
	cout << "E2: " << ques2[1] << endl << endl;
	cout << "Question 3" << endl;
	double* partc[11];
	double* partc_2[11];
	double parta[11];
	double partb[11];
	double prices[1100];
	double prices_2[1100];
	question3(seed, 15,0.5,20,0.04,0.25, prices, prices_2,parta, partb, partc, partc_2);
	string greek[5] = { "Delta ", "Gamma ", "Theta ", "Vega ", "Rho " };
	cout << "monte carlo prices: ";
	for (int i = 0; i < 11; ++i) {
		cout << parta[i] << ", ";
	}
	cout << "\n";
	cout << "black schole prices: ";
	for (int i = 0; i < 11; ++i) {
		cout << partb[i] << ", ";
	}
	cout << "\n";
	
	//ofstream myfile;
	//myfile.open("greeks.csv", fstream::out);
	for (int i = 0; i < 5; ++i) {
		/*cout << greek[i] << "BS ";
		for (int j = 0; j < 11; ++j) {
			cout << partc[j][i] << " ";
		}
		cout << "\n";*/
		cout << greek[i] << "Sim: ";
		for (int j = 0; j < 11; ++j) {
			cout << partc_2[j][i] << ", ";
			//myfile << to_string(partc_2[j][i]) << "," << to_string(partc_2[j][i+1]) << "," << to_string(partc_2[j][i+2]) << "," << to_string(partc_2[j][i+3]) << "," << to_string(partc_2[j][i+4]) << "\n";
		}
		//myfile.close();
		cout << "\n";
	}

	cout << "\n";
	cout << "Question 4" << endl;
	double* ques4 = question4(seed);
	cout << "Reflection " << ques4[0] << endl;
	cout << "Partial " << ques4[1] << endl;
	cout << "Full " << ques4[2] << endl << endl;

	cout << "Question 5" << endl;
	double q5parta[200];
	double q5partb[200];
	double q5partc[200];
	double* ques5 = question5(seed, q5parta, q5partb, q5partc);
	/*ofstream myfile;
	myfile.open("halton.csv", fstream::out);
	for (int i = 0; i < 100; ++i) {
		myfile << to_string(q5parta[i]) << "," << to_string(q5parta[i+100]) << "," << to_string(q5partb[i]) << "," << to_string(q5partb[i + 100]) << "," << to_string(q5partc[i]) << "," << to_string(q5partc[i+ 100]) << "\n";
	}
	myfile.close();*/
	cout << "base (2,7) " << ques5[0] << endl;
	cout << "base (2,4) " << ques5[1] << endl;
	cout << "base (5,7) " << ques5[2] << endl;
}


