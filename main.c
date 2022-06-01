#include "headers.h"

int main(int argc, char **argv) {
    int rank, size, n, m, k;
    char* filename;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    switch (argc)
    {
        case 4: break;
        case 5: filename = argv[4]; break;
        default: return fail;
    }

    n=atoi(argv[1]);
    m=atoi(argv[2]);
    k=atoi(argv[3]);

    if (n < 1 || m > n || (argc == 4 && (k < 0 || k > 4)) || (argc == 5 && k != 0))
    {
        MPI_Finalize();
        return fail;
    }

    double **matrix, **ident_matrix;
    double *matrix_head, *ident_matrix_head;
    double *buf = (double *) malloc(n * sizeof(double));
    double *buf2 = (double *) malloc(n * sizeof(double));
    int *perm = (int *) malloc(n * sizeof(int));
    int *source = (int *) malloc(n * sizeof(int));
    int sizetmp = (n/size) + (rank < n % size);
    int shift = rank * (n/size) + min_(n % size, rank);
    struct timespec start, end;

    allocate_memory(&matrix, &matrix_head, sizetmp, n);
    allocate_memory(&ident_matrix, &ident_matrix_head, sizetmp, n);
    fullfil_source(source, n, size);

    if (argc == 4)
        input(matrix, n, k, rank, source, shift);
    else
        if (file_input(matrix, buf, filename, n, rank, source, shift) == fail)
            goto Bad_Final;

    print(matrix, buf, n, n, m, rank, source, shift);
    if (rank == 0)
        printf("\n");

    timespec_get(&start, TIME_UTC);
    if (solve(matrix, ident_matrix, buf, buf2, perm, source, n, rank, size) == fail)
        goto Bad_Final;
    timespec_get(&end, TIME_UTC);
    double elapsed_time = (double)(((end.tv_sec - start.tv_sec) * 10e9 + (end.tv_nsec - start.tv_nsec)) / 10e9);

    if (argc == 4)
        input(matrix, n, k, rank, source, shift);
    else
        file_input(matrix, buf, filename, n, rank, source, shift);

    if (m == -1)
    {
        print(ident_matrix, buf, n, n, n, rank, source, shift);
        goto Good_Final;
    }

    double res = residual(matrix, ident_matrix, buf, source,  n, rank, shift, sizetmp);

    if (m == 0)
    {
        if (rank == 0)
        {
            printf("residual is %10.3e\ntime elapsed %.2lf sec\n", res, elapsed_time);
        }
        goto Good_Final;
    }

    print(ident_matrix, buf, n, n, m, rank, source, shift);
    if (rank == 0) {
        printf("\nresidual is %10.3e\ntime elapsed %.2lf sec\n", res, elapsed_time);
    }

    Good_Final:
    free(matrix_head);
    free(ident_matrix_head);
    free(matrix);
    free(ident_matrix);
    free(buf);
    free(perm);
    MPI_Finalize();
    return ok;

    Bad_Final:
    free(matrix_head);
    free(ident_matrix_head);
    free(matrix);
    free(ident_matrix);
    free(buf);
    free(perm);
    MPI_Finalize();
    return fail;
}