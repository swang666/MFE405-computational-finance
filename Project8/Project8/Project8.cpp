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
	
	cout << vasicek_model_zero_bond(0.05, 0.18, 0.82, 0.05, 1000, 0.5, 0, 123456) << endl;
	cout << vasicek_model_semicoupon_bond(0.05, 0.18, 30, 0.82, 0.05, 1000, 4, 0, 12345) << endl;
	cout << vasicek_model_zero_bond_call(0.05, 0.18, 0.82, 980, 0.05, 1000, 0.25, 0.5, 0, 12345) << endl;
	cout << cir_model_zero_bond(0.05, 0.18, 0.92, 980, 0.055, 1000, 0.5, 1, 0, 12345) << endl;
	cout << cir_explicit(0.05, 0.18, 0.92, 980, 0.055, 1000, 0.5, 1, 0) << endl;
	cout << G2pp_mc(0.03, 0, 0, 0.03, 0.7, 0.1, 0.3, 0.03, 0.08, 0.03, 985, 1000, 0.5, 1, 0, 12345) << endl;
}
