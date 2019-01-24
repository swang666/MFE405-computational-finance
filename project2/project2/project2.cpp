// project2.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
	cout << question1(seed) << endl << endl;
	cout << "Question 2" << endl;
	cout << question2(seed) << endl << endl;
	cout << "Question 3" << endl;
	cout << "a" << endl;
	cout << "Ea1: " << question3(seed)[0] << endl;
	cout << "Variance is " << question3(seed)[8] << endl;
	cout << "Ea2: " << question3(seed)[1] << endl;
	cout << "Variance is " << question3(seed)[10] << endl;
	cout << "Ea3: " << question3(seed)[2] << endl;
	cout << "Variance is " << question3(seed)[11] << endl;
	cout << "Ea4: " << question3(seed)[3] << endl;
	cout << "Variance is " << question3(seed)[12] << endl;
	cout << "c" << endl;
	cout << "Ec1: " << question3(seed)[4] << endl;
	cout << "Variance is " << question3(seed)[9] << endl;
	cout << "Ec2: " << question3(seed)[5] << endl;
	cout << "Variance is " << question3(seed)[13] << endl;
	cout << "Ec3: " << question3(seed)[6] << endl;
	cout << "Variance is " << question3(seed)[14] << endl;
	cout << "Ec4: " << question3(seed)[7] << endl;
	cout << "Variance is " << question3(seed)[15] << endl << endl;
	cout << "Question 4" << endl;
	cout << "Ca1: " << question4(seed)[0] << endl;
	cout << "Variance is " << question4(seed)[3] << endl;
	cout << "Black Schole: " << question4(seed)[1] << endl;
	cout << "Cb2: " << question4(seed)[2] << endl;
	cout << "Variance is " << question4(seed)[4] << endl << endl;
	cout << "Question 5" << endl;
	double parta[10];
	double* partb[6];
	double partc[10];
	double* partd[6];
	question5(seed, parta, partb, partc, partd);
	cout << "a" << endl;
	for (int i = 0; i < 10; ++i) {
		cout << parta[i] << ", ";
	}
	cout << "\n" << "d" << endl;
	for (int i = 0; i < 10; ++i) {
		cout << partc[i] << ", ";
	}
	/*ofstream myfile;
	myfile.open("ques5a.csv", fstream::out);
	for (int i = 0; i < 10; ++i) {
		myfile << to_string(parta[i]) << "," << to_string(partc[i]) << "\n";
	}
	myfile.close();
	myfile.open("ques5b.csv", fstream::out);
	for (int i = 0; i < 1000; ++i) {
		for (int j = 0; j < 6; ++j) {
			myfile << to_string(partb[j][i]) << ",";
		}
		myfile << "\n";
	}
	myfile.close();
	myfile.open("ques5d.csv", fstream::out);
	for (int i = 0; i < 1000; ++i) {
		for (int j = 0; j < 6; ++j) {
			myfile << to_string(partd[j][i]) << ",";
		}
		myfile << "\n";
	}
	myfile.close();*/
	cout << "\n" << endl;
	cout << "Question 6" << endl;
	cout << "6a: " << question6(seed)[0] << endl;
	cout << "6b: " << question6(seed)[1] << endl;
	cout << "Variance is " << question6(seed)[3] << endl;
	cout << "6c: " << question6(seed)[2] << endl;
	cout << "Variance is " << question6(seed)[4] << endl;
}

