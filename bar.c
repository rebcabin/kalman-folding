#include <stdio.h>

void print_RMO_matrix (double * M, int n_rows, int n_cols) {
    for (int row = 0; row < n_rows; ++row)
    {   for (int col = 0; col < n_cols; ++col)
        {   printf ("%g", M[col + n_cols * row]);
            if (n_cols-1 != col)
            {   printf (", ");   }   }
        if (n_rows-1 != row)
        {   printf ("\n");   }   }
    printf ("\n\n");   }

void print_CMO_matrix (double * M, int n_rows, int n_cols) {
    for (int row = 0; row < n_rows; ++row)
    {   for (int col = 0; col < n_cols; ++col)
        {   printf ("%g", M[row + n_rows * col]);
            if (n_cols-1 != col)
            {   printf (", ");   }   }
        if (n_rows-1 != row)
        {   printf ("\n");   }   }
    printf ("\n\n");   }

int main (int argc, char ** argv) {
    int N = 3;
    int M = 3;
    int MN = 3 * 3;
    /* The following input appears in row-major order. The underlying fortran
       code 'sees' the transpose of this matrix, i.e., the same data in
       column-major order. Although this is the matrix A for fortran from our
       example, it looks, in C, like AT.

       BLAS does not include a high-performance matrix transpose routine. */

    double A[4][5] = { { 1.,  1.,  1.,  1.,  1.},
                       { 0.,  1., -1., -2.,  2.},
                       { 0.,  1.,  1.,  4.,  4.},
                       { 0.,  1., -1., -8.,  8.}};

    double DX[3][3] = { {1 , 2 ,  3},
                        {4 , 5 ,  6},
                        {7 , 8 , 19} };
    double DY[3][3] = { {0 , 0 ,  0},
                        {0 , 0 ,  0},
                        {0 , 0 ,  0} };

    print_RMO_matrix ((double *)A, 4, 5);
    print_CMO_matrix ((double *)A, 5, 4);

    int INCX = 1;
    int INCY = 1;
    /* from http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f_source.html
       SUBROUTINE dcopy(N,DX,INCX,DY,INCY)
       *
       *  -- Reference BLAS level1 routine (version 3.4.0) --
       *  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
       *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
       *     November 2011
       *
       *     .. Scalar Arguments ..
       INTEGER INCX,INCY,N
       *     ..
       *     .. Array Arguments ..
       DOUBLE PRECISION DX(*),DY(*)
       *     ..
       */

    extern void dcopy_ (int * N, double * DX, int * INCX, double * DY, int * INCY);

    print_RMO_matrix ((double *)DX, N, N);
    print_RMO_matrix ((double *)DY, N, N);

    dcopy_ (&MN, DX[0], &INCX, DY[0], &INCY);

    print_RMO_matrix ((double *)DX, N, N);
    print_RMO_matrix ((double *)DY, N, N);   }
