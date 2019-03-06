// Project8.cpp : This file contains the 'main' function. Program execution begins and ends there.
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

	//cout << vasicek_model_zero_bond(0.05, 0.18, 0.82, 0.05, 1000, 0.5, 0, 123456) << endl;
	cout << cir_model_zero_bond(0.05, 0.18, 0.92, 980, 0.055, 1000, 0.5, 1, 0, 12345) << endl;
	cir_explicit(0.05, 0.18, 0.92, 980, 0.055, 1000, 0.5, 1, 0);
}
