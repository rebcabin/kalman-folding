#include <stdio.h>
#include <gsl/gsl_blas.h>





int main() {
#include <stddef.h>
#include <lapacke.h>
int N = 3;
int NN = 9;
double M[3][3] = { {1 , 2 , 3},
                   {4 , 5 , 6},
                   {7 , 8 , 9} };
int pivotArray[3]; //since our matrix has three rows
int errorHandler;
double lapackWorkspace[9];

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
