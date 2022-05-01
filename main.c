#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "itermethods.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: \n\t-a for Jacobi method. \n\t-o for Jordan method. \n\t-r for relaxation method.\n");
        exit(1);
    }

    FILE *DATA, *RESULT;
    DATA = fopen("data.dat", "r");
    RESULT = fopen("result.dat", "w+");

    int n;
    if (getc(DATA) != '#')
    {
        fprintf(stderr, "Error: file in wrong format.");
        exit(1);
    }

    fscanf(DATA, "%d", &n);

    double **A = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i)
    {
        A[i] = (double *)malloc(sizeof(double) * (n + 1));
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            fscanf(DATA, "%lf", &A[i][j]);
        }
    }
    for (int i = 0; i < n; ++i)
    {
        fscanf(DATA, "%lf", &A[i][n]);
    }

    double *X = (double *) malloc(sizeof(double) * n);

    int opt;
    while ((opt = getopt(argc, argv, "aor")) != -1)
    {
        switch (opt)
        {
        case 'a':
            solve_jacoby(A, X);
            break;
        case 'o':
            solve_jordan(A, X);
            break;
        case 'r':
            solve_relaxation(A, X);
            break;
        default:
            fprintf(stderr, "Usage: \n\t-a for Jacobi method. \n\t-o for Jordan method. \n\t-r for relaxation method.\n");
            exit(1);
        }
        break;
    }

    fprintf(RESULT, "# %d \n", n);

    for (int i = 0; i < n; ++i)
    {
        fprintf(RESULT, "%10e\n", X[i]);
    }

    return 0;
}