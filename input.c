#include "headers.h"

void input(double** x, int n, int k, int rank, int *source, int shift)
{
    for (int i=0; i<n; ++i)
    {
        if (source[i] == rank)
            for (int j=0; j<n; ++j)
            {
                x[i-shift][j]=f(n, k, i+1, j+1);
            }
    }
}

double f(int n, int k, int i, int j)
{
    switch (k)
    {
        case 1:return (double)(n-max_(i,j)+1);
        case 2:return (double)(max_(i,j));
        case 3:return (double)abs(i-j);
        case 4:return 1.0/(double)(i+j-1);
        default: return 0.0;
    }
}

int file_input(double** x, double *buf, const char* filename, int n, int rank, int* source, int shift)
{
    FILE *f;
    if (rank == 0){
        f=fopen(filename,"r");
        if (!f)
        {
            return fail;
        }
    }
    for (int i=0; i<n; ++i)
    {
        if (rank == 0){
            for (int j=0; j<n; ++j)
            {
                if (fscanf(f,"%lf",&buf[j]) != 1)
                {
                    fclose(f);
                    return fail;
                }
            }
            if (source[i] > 0)
                MPI_Send(buf, n, MPI_DOUBLE, source[i], 0, MPI_COMM_WORLD);
        }

        if (rank == source[i])
        {
            if (rank > 0)
                MPI_Recv(buf, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            copy_raw(x[i - shift], buf, n);
        }
    }
    if (rank == 0)
        fclose(f);
    return ok;
}
