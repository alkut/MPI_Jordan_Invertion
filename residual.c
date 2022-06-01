#include "headers.h"

double residual(double **matrix, double **ident_matrix, double *buf, int *source, int n, int rank, int shift, int sizetmp)
{
    double ans = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double tmp = -(double)(i == j);
            tmp += scalar_prod(matrix, ident_matrix, buf, n, i, j, rank, source, shift, sizetmp);
            ans += tmp * tmp;
        }
    }
    return sqrt(ans);
}