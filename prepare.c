#include "headers.h"

void allocate_memory(double ***matrix, double **head, int sizetmp, int n)
{
    *matrix = (double **)malloc(sizetmp * sizeof(double *));
    double *mem = (double *) malloc(sizetmp * n * sizeof(double));
    *head = mem;
    for (int i=0; i < sizetmp; ++i)
    {
        (*matrix)[i] = mem + i*n;
    }
}

void fullfil_source(int *source, int n, int size)
{
    int bsize = (n / size) + 1;
    int ssize = n / size;
    for (int i = 0; i < n % size; ++i)
    {
        for (int j = 0; j < bsize; ++j)
        {
            source[i * bsize + j] = i;
        }
    }
    for (int i = n % size; i < size; ++i)
    {
        for (int j = 0; j < ssize; ++j)
        {
            source[bsize * (n % size) + (i - n % size) * ssize + j] = i;
        }
    }
}

void fullfil_perm(int *perm, int n)
{
    for (int i = 0; i < n; ++i)
    {
        perm[i] = i;
    }
}

void make_ident_matrix(double **matr, int sizetmp, int n, int shift)
{
    for (int i = 0; i < sizetmp; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            matr[i][j] = (j == i + shift);
        }
    }
}