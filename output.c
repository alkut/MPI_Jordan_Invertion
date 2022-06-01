#include "headers.h"

void print(double** x, double *buf, int m, int n, int k, int rank, int *source, int shift)
{
    for (int i=0; i<min_(m, k); ++i)
    {
        if (source[i] == rank)
        {
            copy_raw(buf, x[i-shift], min_(n, k));
            if (rank > 0)
                MPI_Send(buf, min_(n, k), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 0) {
            if (source[i] > 0)
                MPI_Recv(buf, min_(n, k), MPI_DOUBLE, source[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < min_(n, k); ++j) {
                printf("%10.3e ", buf[j]);
            }
            printf("\n");
        }
    }
}
