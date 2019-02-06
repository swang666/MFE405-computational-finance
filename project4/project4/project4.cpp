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
	
}
