#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define fail  -1
#define ok  0
#define machine_eps  2.20e-16

#define max_(i,j) (i<j ? j: i)
#define min_(i,j) (i<j ? i: j)

typedef struct dint
{
    double d;
    int i;
} dint;

//--------begin-input-----------------------------------------------------------------------------------------//
void input(double** x, int n, int k, int rank, int *source, int shift);
double f(int n, int k, int i, int j);
int file_input(double** x, double *buf, const char* filename, int n, int rank, int* source, int shift);
//--------end---input-----------------------------------------------------------------------------------------//

//--------begin-output----------------------------------------------------------------------------------------//
void print(double** x, double *buf, int m, int n, int k, int rank, int *source, int shift);
//--------end---output----------------------------------------------------------------------------------------//

//--------begin-utilities-------------------------------------------------------------------------------------//
void swap_raw(double **matrix, double *buf, int n, int rank, int raw1, int raw2, int *source, int shift);
void swap_columns(double **matrix, int *source, int n, int rank, int col1, int col2, int shift);
void copy_col(double **matr_to, double **matr_from, int to, int from, int n, int rank, int *source, int shift);
double scalar_prod(double **matr1, double **matr2, double *buf, int n, int raw1, int col2, int rank, int *source, int shift, int sizetmp);
void copy_raws(double **matr_to, double **matr_from, double *buf, int to, int from, int n, int rank, int *source, int shift);
void copy_raw(double *raw_to, double *raw_from, int n);
void swap(double *x, double *y);
void swapp(int* x, int* y);
void swappp(double** x, double **y);
void swapppp(double *x, double *y, int n);
//--------end---utilities-------------------------------------------------------------------------------------//

//--------begin-jordan-invertion------------------------------------------------------------------------------//
int solve(double **matrix, double **ident_matrix, double *buf, double *buf2, int *perm, int *source, int n, int rank, int size);
int jordan_invertion(double **matrix, double **ident_matrix, double *buf, double *buf2, int *perm, int sizetmp, int n, int step, int rank, int *source, int shift);
int make_zero(double **matrix, double **ident_matrix, double *buf, double *buf2, int *perm, int sizetmp, int n, int step, int rank, int *source, int shift);
void move_max(double **matrix, double **ident_matrix, double *buf, int *perm, int sizetmp, int n, int step, int rank, int *source, int shift);
void get_max(double **matrix, int *source, int n, int rank, int step, int shift, dint *loc);
//--------end---jordan-invertion------------------------------------------------------------------------------//

//--------begin-residual--------------------------------------------------------------------------------------//
double residual(double **matrix, double **ident_matrix, double *buf, int *source, int n, int rank, int shift, int sizetmp);
//--------end---residual--------------------------------------------------------------------------------------//

//--------begin-prepare---------------------------------------------------------------------------------------//
void allocate_memory(double ***matrix, double **head, int sizetmp, int n);
void fullfil_source(int *source, int n, int size);
void fullfil_perm(int *perm, int n);
void make_ident_matrix(double **matr, int sizetmp, int n, int shift);
//--------end---prepare---------------------------------------------------------------------------------------//

