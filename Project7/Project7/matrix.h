#ifndef MATRIX_H
#define MATRIX_H
#include <stdint.h>

double** zero_matrix(int n, int m);
void print_matrix(double** mat, int n, int m);
double** matrix_mult(double** mat1, double** mat2, int n1, int m1, int n2, int m2);
double** matrix_add(double** mat1, double** mat2, int n1, int m1);
double** back_sub(double** A, double** b, int n);
double** column_bind(double** c1, double** c2, int n1, int m1);
double** upper_tri_inverse(double** A, int n);
#endif //MATRIX_H