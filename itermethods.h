#ifndef ITERMETHODS_H
#define ITERMETHODS_H

int solve_jacoby(double **A, double *X, int n);
int solve_jordan(double **A, double *X, int n);
int solve_relaxation(double **A, double *X, int n);
int is_diagonally_dominant(double **A, int n);
double get_resedue(double *X1, double *X2, int n);

#endif // ITERMETHODS_H