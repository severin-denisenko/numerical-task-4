#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "itermethods.h"

/*
    Вычисление прерывается если ||X1 - X2|| < min_tolerance или k > max_iterations.
*/

double min_tolerance = 10e-10;
int max_iterations = 10000;

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

    for (int i = 0; i < n; i++)
    {
        X1[i] = 0;
        X2[i] = 0;
    }

    for (int k = 0; k < max_iterations; k++)
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

        // Debug
        // for (int i = 0; i < n; i++)
        // {
        //     printf("%le ", X2[i]);
        // }
        // printf("\n");

        if (get_residue(X1, X2, n) < min_tolerance)
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

    for (int i = 0; i < n; i++)
    {
        X1[i] = 0;
        X2[i] = 0;
    }

    for (int k = 0; k < max_iterations; k++)
    {
        for (int i = 0; i < n; i++)
        {
            X2[i] = A1[i][n];
            for (int j = i + 1; j < n; j++)
            {
                X2[i] -= A1[i][j] * X1[j];
            }

            for (int j = 0; j < i; j++)
            {
                X2[i] -= A1[i][j] * X2[j];
            }
        }

        if (get_residue(X1, X2, n) < min_tolerance)
        {
            break;
        }

        // Debug
        // for (int i = 0; i < n; i++)
        // {
        //     printf("%le ", X2[i]);
        // }
        // printf("\n");

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
    // Копирование с делением

    double **P = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i)
    {
        P[i] = (double *)malloc(sizeof(double) * n);
    }

    double *Q1 = (double *)malloc(sizeof(double) * n);
    double *Q2 = (double *)malloc(sizeof(double) * n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (A[i][i] == 0)
            {
                fprintf(stderr, "Error: element (%d, %d) is 0.\n", i, i);
                exit(1);
            }
            P[i][j] = - A[i][j] / A[i][i];
        }
    }

    for (int i = 0; i < n; i++)
    {
        Q1[i] = A[i][n] / A[i][i];
        Q2[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        X[i] = 0;
    }

    for (int k = 0; k < max_iterations; k++)
    {
        int max_i = 0;

        for (int i = 0; i < n; i++)
        {
            if (fabs(Q1[max_i]) < fabs(Q1[i]))
            {
                max_i = i;
            }
        }

        if (fabs(Q1[max_i]) < min_tolerance)
        {
            break;
        }

        X[max_i] = X[max_i] + Q1[max_i];

        for (int i = 0; i < n; i++)
        {
            Q2[i] = Q1[i] + P[i][max_i] * Q1[max_i];
        }

        // Debug
        // for (int i = 0; i < n; i++)
        // {
        //     printf("%le ", X[i]);
        // }
        // printf("\n");

        double *tmp = Q1;
        Q1 = Q2;
        Q2 = tmp;
    }

    for (int i = 0; i < n; i++)
    {
        free(P[i]);
    }
    free(P);
    free(Q1);
    free(Q2);

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

double get_residue(double *X1, double *X2, int n)
{
    double r = 0;

    for (int i = 0; i < n; i++)
    {
        r += (X1[i] - X2[i]) * (X1[i] - X2[i]);
    }

    return sqrt(r);
}
