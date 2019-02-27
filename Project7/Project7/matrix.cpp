#include <iostream>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>  
#include "matrix.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace std;

const double pi = atan(1) * 4;

double** zero_matrix(int n, int m) {
	double** out = new double*[n];
	for (int i = 0; i < n; ++i) {
		out[i] = new double[m];
		for (int j = 0; j < m; ++j) {
			out[i][j] = 0;
		}
	}
	return out;
}

void print_matrix(double** mat, int n, int m) {
	for (int i = 0; i < n;++i) {
		for (int j = 0; j < m; ++j) {
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
}

double** matrix_mult(double** mat1, double** mat2, int n1, int m1, int n2, int m2) {
	//m1 should be equal to n2
	double** out = new double*[n1];
	for (int i = 0; i < n1; ++i) {
		out[i] = new double[m2];
		for (int j = 0; j < m2; ++j) {
			double sum = 0;
			for (int k = 0; k < m1; ++k) {
				sum = sum + mat1[i][k] * mat2[k][j];
			}
			out[i][j] = sum;
		}
	}
	return out;
	
}

double** matrix_add(double** mat1, double** mat2, int n1, int m1) {
	// n1 = n2, m1 = m2
	double** out = new double*[n1];
	for (int i = 0; i < n1; ++i) {
		out[i] = new double[m1];
		for (int j = 0; j < m1; ++j) {
			out[i][j] = mat1[i][j] + mat2[i][j];
		}
	}
	return out;
}

double** back_sub(double** A, double** b, int n) {
	double sum;
	double** x = zero_matrix(n, 1);
	for (int i = n-1; i >= 0; --i) {
		sum = 0;
		for (int j = i; j < n - 1; ++j) {
			sum = sum + A[i][j+1] * x[j+1][0];
		}
		x[i][0] = (b[i][0] - sum) / A[i][i];
	}
	return x;
}

double** column_bind(double** c1, double** c2, int n1, int m1) {
	double** out = zero_matrix(n1, m1 + 1);
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < m1; ++j) {
			out[i][j] = c1[i][j];
		}
		out[i][m1] = c2[i][0];
	}
	return out;
}

double** upper_tri_inverse(double** A, int n) {
	double** e;
	double** x;
	double** temp;
	e = zero_matrix(n, 1);
	e[0][0] = 1;
	x = back_sub(A, e, n);
	for (int i = 1; i < n; ++i) {
		e = zero_matrix(n, 1);
		e[i][0] = 1;
		temp = back_sub(A, e, n);
		x = column_bind(x, temp, n, i);
	}
	return x;
}