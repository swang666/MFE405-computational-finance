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
	int64_t seed = 1234567;
	cout << "Question 1" << endl;

	double * ques1 = question1(seed);
	cout << "Prob: " << ques1[0] << endl;
	cout << "E1: " << ques1[1] << endl;
	cout << "E2: " << ques1[2] << endl;
	cout << "E3: " << ques1[3] << endl << endl;
	cout << "Question 2" << endl;
	double* ques2 = question2(seed);
	cout << ques2[0] << endl;
	cout << ques2[1] << endl << endl;
	cout << "Question 3" << endl;
	double* partc[11];
	double* partc_2[11];
	double parta[2];
	question3(seed, 15,0.5,20,0.04,0.25,parta, partc, partc_2);
	string greek[5] = { "Delta ", "Gamma ", "Theta ", "Vega ", "Rho " };
	cout << parta[0] << endl;
	cout << parta[1] << endl;
	for (int i = 0; i < 5; ++i) {
		cout << greek[i] << "BS ";
		for (int j = 0; j < 11; ++j) {
			cout << partc[j][i] << " ";
		}
		cout << "\n";
		cout << greek[i] << "Sim ";
		for (int j = 0; j < 11; ++j) {
			cout << partc_2[j][i] << " ";
		}
		cout << "\n";
	}


}


