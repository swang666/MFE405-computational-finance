// project4.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "binom.h"
#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;

int main()
{
	double** ques1 = question1();
	cout << "1a" << endl;
	for (int i = 0; i < 7; ++i) {
		cout << ques1[0][i] << endl;
	}
	cout << "1b" << endl;
	for (int i = 0; i < 7; ++i) {
		cout << ques1[1][i] << endl;
	}
	cout << "1c" << endl;
	for (int i = 0; i < 7; ++i) {
		cout << ques1[2][i] << endl;
	}
	cout << "1d" << endl;
	for (int i = 0; i < 7; ++i) {
		cout << ques1[3][i] << endl;
	}
	cout << endl;
	cout << "Question 2" << endl;
	cout << question2() << endl << endl;
	cout << "Question 3" << endl;
	double** ques3 = question3();
	for (int i = 0; i < 31; ++i) {
		cout << ques3[0][i] << " ";
	}
	cout << endl << endl;
	for (int i = 0; i < 31; ++i) {
		cout << ques3[1][i] << " ";
	}
	cout << endl << endl;
	for (int i = 0; i < 31; ++i) {
		cout << ques3[2][i] << " ";
	}
	cout << endl << endl;
	for (int i = 0; i < 31; ++i) {
		cout << ques3[3][i] << " ";
	}
	cout << endl << endl;
	for (int i = 0; i < 31; ++i) {
		cout << ques3[4][i] << " ";
	}
	cout << endl << endl;
	for (int i = 0; i < 38; ++i) {
		cout << ques3[5][i] << " ";
	}
	cout << endl << endl;
	cout << "Question 4" << endl;
	double** ques4 = question4();
	for (int i = 0; i < 11; ++i) {
		cout << ques4[0][i] << " ";
	}
	cout << endl << endl;
	for (int i = 0; i < 11; ++i) {
		cout << ques4[1][i] << " ";
	}
	cout << endl << endl;
	cout << "Question 5" << endl;
	double** ques5 = question5();
	cout << "5a" << endl;
	for (int i = 0; i < 9; ++i) {
		cout << ques5[0][i] << endl;
	}
	cout << "5b" << endl;
	for (int i = 0; i < 9; ++i) {
		cout << ques5[1][i] << endl;
	}
}
