#include "headers.h"

void swap_raw(double **matrix, double *buf, int n, int rank, int raw1, int raw2, int *source, int shift)
{
    if (raw1 == raw2) return;
    if (source[raw1] == source[raw2])
    {
        if (rank == source[raw1])
            swappp(&matrix[raw1 - shift], &matrix[raw2 - shift]);
        return;
    }
    if (rank == source[raw1])
    {
        copy_raw(buf, matrix[raw1 - shift], n);
        MPI_Send(buf, n, MPI_DOUBLE, source[raw2], 0, MPI_COMM_WORLD);
        MPI_Recv(buf, n, MPI_DOUBLE, source[raw2], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        copy_raw(matrix[raw1 - shift], buf, n);
    }
    if (rank == source[raw2])
    {
        MPI_Recv(buf, n, MPI_DOUBLE, source[raw1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        swapppp(buf, matrix[raw2 - shift], n);
        MPI_Send(buf, n, MPI_DOUBLE, source[raw1], 0, MPI_COMM_WORLD);
    }
}

void swap_columns(double **matrix, int *source, int n, int rank, int col1, int col2, int shift)
{
    if (col1 == col2) return;
    for (int i = 0; i < n; ++i)
    {
        if (source[i] == rank)
        {
            swap(matrix[i - shift] + col1, matrix[i - shift] + col2);
        }
    }
}

void copy_col(double **matr_to, double **matr_from, int to, int from, int n, int rank, int *source, int shift)
{
    for (int i = 0; i < n; ++i)
    {
        if (rank == source[i])
        {
            matr_to[i - shift][to] = matr_from[i - shift][from];
        }
    }
}

double scalar_prod(double **matr1, double **matr2, double *buf, int n, int raw1, int col2, int rank, int *source, int shift, int sizetmp)
{
    if (rank == source[raw1])
        copy_raw(buf, matr1[raw1 - shift], n);
    MPI_Bcast(buf, n, MPI_DOUBLE, source[raw1], MPI_COMM_WORLD);
    double loc = 0.0, gl;
    for (int i = 0; i < sizetmp; ++i)
    {
        loc += buf[i + shift] * matr2[i][col2];
    }
    MPI_Allreduce(&loc, &gl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return gl;
}

void copy_raws(double **matr_to, double **matr_from, double *buf, int to, int from, int n, int rank, int *source, int shift)
{
    if (source[to] == source[from])
    {
        if (rank == source[to])
        {
            copy_raw(matr_to[to - shift], matr_from[from - shift], n);
        }
        return;
    }

    if (rank == source[from])
    {
        copy_raw(buf, matr_from[from - shift], n);
        MPI_Send(buf, n, MPI_DOUBLE, source[to], 0, MPI_COMM_WORLD);
    }
    if (rank == source[to])
    {
        MPI_Recv(buf, n, MPI_DOUBLE, source[from], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        copy_raw(matr_to[to - shift], buf, n);
    }
}

void copy_raw(double *raw_to, double *raw_from, int n)
{
    for (int j=0; j<n; ++j)
    {
        raw_to[j] = raw_from[j];
    }
}

void swap(double *x, double *y)
{
    double tmp=*x;
    *x=*y;
    *y=tmp;
}

void swapp(int* x, int* y)
{
    unsigned tmp=*x;
    *x=*y;
    *y=tmp;
}

void swappp(double **x, double **y)
{
    double *tmp= *x;
    *x=*y;
    *y=tmp;
}

void swapppp(double *x, double *y, int n)
{
    for (int i = 0; i< n; ++i)
    {
        swap(x+i, y+i);
    }
}
