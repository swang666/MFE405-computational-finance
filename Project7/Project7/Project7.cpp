// Project7.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "matrix.h"
#include "fdm.h"
#include <iostream>

using namespace std;

int main()
{
	cout << EFDM(10, 1) << endl;
	cout << EFDM(10, 3) << endl;
	cout << EFDM(10, 4) << endl;
	double** A = zero_matrix(3, 3);
	double** b = zero_matrix(3, 1);
	A[0][0] = 1;
	A[0][1] = 2;
	A[0][2] = 4;
	A[1][1] = 5;
	A[1][2] = 8;
	A[2][2] = 9;
	double** Ainv = upper_tri_inverse(A, 3);
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			cout << Ainv[i][j] << " ";
		}
		cout << endl;
	}
}
