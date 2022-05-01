#ifndef ITERMETHODS_H
#define ITERMETHODS_H

int solve_jacoby(double **A, double *X);
int solve_jordan(double **A, double *X);
int solve_relaxation(double **A, double *X);
int check_for_diagonally_dominant(double **A);

#endif // ITERMETHODS_H