#include "headers.h"

int solve(double **matrix, double **ident_matrix, double *buf, double *buf2, int *perm, int *source, int n, int rank, int size)
{
    int sizetmp = (n/size) + (rank < n % size);
    int shift = rank * (n/size) + min_(n % size, rank);

    fullfil_perm(perm, n);

    make_ident_matrix(ident_matrix, sizetmp, n, shift);

    return jordan_invertion(matrix, ident_matrix, buf, buf2, perm, sizetmp, n, 0, rank, source, shift);

}

int jordan_invertion(double **matrix, double **ident_matrix, double *buf, double *buf2, int *perm, int sizetmp, int n, int step, int rank, int *source, int shift)
{
    for (int i = 0; i < n; ++i)
    {
        if (make_zero(matrix, ident_matrix, buf, buf2, perm, sizetmp, n, i, rank, source, shift) == fail)
            return fail;
    }
    for (int i = 0; i < n; ++i)
    {
        copy_col(matrix, ident_matrix, perm[i], i, n, rank, source, shift);
    }
    for (int i = 0; i < n; ++i)
    {
        copy_raws(ident_matrix, matrix, buf, perm[i], i, n, rank, source, shift);
    }
    return ok;
}

int make_zero(double **matrix, double **ident_matrix, double *buf, double *buf2, int *perm, int sizetmp, int n, int step, int rank, int *source, int shift)
{
    move_max(matrix, ident_matrix, buf, perm, sizetmp, n, step, rank, source, shift);
    if (rank == source[step])
    {
        double coef = matrix[step - shift][step];
        if (fabs(coef) < machine_eps)
            return fail;
        for (int i = 0; i < n; ++i)
        {
            matrix[step - shift][i] /= coef;
            ident_matrix[step - shift][i] /= coef;
        }
        copy_raw(buf, matrix[step - shift], n);
        copy_raw(buf2, ident_matrix[step - shift], n);
    }

    MPI_Bcast(buf, n, MPI_DOUBLE, source[step], MPI_COMM_WORLD);
    MPI_Bcast(buf2, n, MPI_DOUBLE, source[step], MPI_COMM_WORLD);

    for (int i = 0; i < n; ++i)
    {
        if (i == step) continue;
        if (rank == source[i])
        {
            double coef = matrix[i - shift][step];
            for (int j = 0; j < n; ++j)
            {
                matrix[i - shift][j] -= buf[j] * coef;
                ident_matrix[i - shift][j] -= buf2[j] * coef;
            }
        }
    }
    return ok;
}

void move_max(double **matrix, double **ident_matrix, double *buf, int *perm, int sizetmp, int n, int step, int rank, int *source, int shift)
{
    dint gl, loc;
    get_max(matrix, source, n, rank, step, shift, &loc);

    MPI_Allreduce(&loc, &gl, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    swap_raw(matrix, buf, n, rank, step, gl.i / n, source, shift);
    swap_raw(ident_matrix, buf, n, rank, step, gl.i / n, source, shift);

    swap_columns(matrix, source, n, rank, step, gl.i % n, shift);
    swap_columns(ident_matrix, source, n, rank, step, gl.i % n, shift);

    swapp(perm + (gl.i % n), perm + step);
}

void get_max(double **matrix, int *source, int n, int rank, int step, int shift, dint *loc)
{
    loc->d = 0.0;
    loc->i = 0;
    for (int i = step; i < n ; ++i)
    {
        if (rank != source[i]) continue;
        for (int j=step; j<n; ++j)
        {
            if (fabs(matrix[i - shift][j]) > loc->d)
            {
                loc->d = fabs(matrix[i - shift][j]);
                loc->i = i * n + j;
            }
        }
    }
}
