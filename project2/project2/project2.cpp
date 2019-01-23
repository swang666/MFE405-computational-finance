// project2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "rng.h"
#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  

using namespace std;

int main()
{
	int64_t seed = 1234567;
	cout << "Question 1" << endl;
	cout << question1(seed) << endl << endl;
	cout << "Question 2" << endl;
	cout << question2(seed) << endl << endl;
	cout << "Question 3" << endl;
	cout << "a" << endl;
	cout << "Ea1: " << question3(seed)[0] << endl;
	cout << "Ea2: " << question3(seed)[1] << endl;
	cout << "Ea3: " << question3(seed)[2] << endl;
	cout << "Ea4: " << question3(seed)[3] << endl;
	cout << "c" << endl;
	cout << "Ec1: " << question3(seed)[4] << endl;
	cout << "Ec2: " << question3(seed)[5] << endl;
	cout << "Ec3: " << question3(seed)[6] << endl;
	cout << "Ec4: " << question3(seed)[7] << endl << endl;
	cout << "Question 4" << endl;
	cout << "Ca1: " << question4(seed)[0] << endl;
}

