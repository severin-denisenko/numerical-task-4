#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "itermethods.h"

/*
    Вычисление перывается если ||X1 - X2|| < TOLERANSE или k > ITERATIONS.
*/

#define TOLERANSE 10e-10
#define ITERATIONS 10000

int solve_jacoby(double **A, double *X, int n)
{
    // Копирование с делением

    double **A1 = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i)
    {
        A1[i] = (double *)malloc(sizeof(double) * (n + 1));
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            if (A[i][i] == 0)
            {
                fprintf(stderr, "Error: element (%d, %d) is 0.\n", i, i);
                exit(1);
            }
            A1[i][j] = A[i][j] / A[i][i];
        }
    }

    double *X1 = (double *)malloc(sizeof(double) * n);
    double *X2 = (double *)malloc(sizeof(double) * n);

    for (int k = 0; k < ITERATIONS; k++)
    {
        for (int i = 0; i < n; i++)
        {
            X2[i] = A1[i][n];
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    continue;
                }

                X2[i] -= A1[i][j] * X1[j];
            }
        }

        if (get_resedue(X1, X2, n) < TOLERANSE)
        {
            break;
        }

        double *tmp = X1;
        X1 = X2;
        X2 = tmp;
    }

    for (int i = 0; i < n; i++)
    {
        X[i] = X2[i];
    }

    free(X1);
    free(X2);

    for (int i = 0; i < n; i++)
    {
        free(A1[i]);
    }
    free(A1);

    return 0;
}

int solve_jordan(double **A, double *X, int n)
{
    // Копирование с делением

    double **A1 = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i)
    {
        A1[i] = (double *)malloc(sizeof(double) * (n + 1));
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            if (A[i][i] == 0)
            {
                fprintf(stderr, "Error: element (%d, %d) is 0.\n", i, i);
                exit(1);
            }
            A1[i][j] = A[i][j] / A[i][i];
        }
    }

    double *X1 = (double *)malloc(sizeof(double) * n);
    double *X2 = (double *)malloc(sizeof(double) * n);

    for (int k = 0; k < ITERATIONS; k++)
    {
        for (int i = 0; i < n; i++)
        {
            X2[i] = A1[i][n];
            for (int j = i + 1; j < n; j++)
            {
                X2[i] -= A1[i][j] * X1[j];
            }

            for (int j = 0; j < i - 1; j++)
            {
                X2[i] -= A1[i][j] * X2[j];
            }
        }
        
        if (get_resedue(X1, X2, n) < TOLERANSE)
        {
            break;
        }

        double *tmp = X1;
        X1 = X2;
        X2 = tmp;
    }

    for (int i = 0; i < n; i++)
    {
        X[i] = X2[i];
    }

    free(X1);
    free(X2);

    for (int i = 0; i < n; i++)
    {
        free(A1[i]);
    }
    free(A1);

    return 0;
}

int solve_relaxation(double **A, double *X, int n)
{
    return 0;
}

int is_diagonally_dominant(double **A, int n)
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                sum += A[j][i];
            }
        }

        if (sum >= A[i][i])
        {
            fprintf(stderr, "Matrix is not diagonally dominant!\n");
            return 0;
        }
    }

    return 1;
}

double get_resedue(double *X1, double *X2, int n)
{
    double r = 0;

    for (int i = 0; i < n; i++)
    {
        r += (X1[i] - X2[i]) * (X1[i] - X2[i]);
    }

    return sqrt(r);
}
