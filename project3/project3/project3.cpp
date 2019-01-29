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
	cout << "E3: " << ques1[3] << endl;


}


