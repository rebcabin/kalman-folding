#include <stdio.h>
#include <stddef.h>





int main() {
int N = 3;
int NN = 9;
double M[3][3] = { {1 , 2 ,  3},
                   {4 , 5 ,  6},
                   {7 , 8 , 19} };
int pivotArray[3]; //since our matrix has three rows
int errorHandler;
double lapackWorkspace[9];
/*
  SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
  *
  *  -- LAPACK routine (version 3.1) --
  *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  *     November 2006
  *
  *     .. Scalar Arguments ..
  INTEGER            INFO, LDA, M, N
  *     ..
  *     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   A( LDA, * )
  */

extern void dgetrf_ (int * m, int * n, double * A, int * LDA, int * IPIV,
                     int * INFO);

/* from http://www.netlib.no/netlib/lapack/double/dgetri.f
  SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
  *
  *  -- LAPACK routine (version 3.1) --
  *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  *     November 2006
  *
  *     .. Scalar Arguments ..
  INTEGER            INFO, LDA, LWORK, N
  *     ..
  *     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   A( LDA, * ), WORK( * )
  */

extern void dgetri_ (int * n, double * A, int * LDA, int * IPIV,
                     double * WORK, int * LWORK, int * INFO);

// dgetrf(M,N,A,LDA,IPIV,INFO) means invert LDA columns of an M by N matrix
// called A, sending the pivot indices to IPIV, and spitting error information
// to INFO. also don't forget (like I did) that when you pass a two-dimensional
// array to a function you need to specify the number of "rows"
dgetrf_(&N, &N, M[0], &N, pivotArray, &errorHandler);
printf ("dgetrf eh, %d, should be zero\n", errorHandler);

dgetri_(&N, M[0], &N, pivotArray, lapackWorkspace, &NN, &errorHandler);
printf ("dgetri eh, %d, should be zero\n", errorHandler);

for (size_t row = 0; row < N; ++row)
 {   for (size_t col = 0; col < N; ++col)
     {   printf ("%g", M[row][col]);
         if (N-1 != col)
         {   printf (", ");   }   }
     if (N-1 != row)
     {   printf ("\n");   }   }
// ~~> produces
return 0;
}
