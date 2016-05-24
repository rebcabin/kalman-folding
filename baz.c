#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>
int main (int argc, char ** argv)
{ const int    b = 1;
  const int    n = 4;
  double x[n * 1] = {0., 0., 0., 0};
  double Z[b * b] = {1.};
  double P[n * n] = {1000.,    0.,    0.,    0.,
                        0., 1000.,    0.,    0.,
                        0.,    0., 1000.,    0.,
                        0.,    0.,    0., 1000. };

  double A[b * n] = {1., 0., 0., 0};
  double z[b] = {-2.28442};

  /* Transcribe the following Wolfram code (the intermediate matrices are not
   * necessary in Wolfram, but we need them in C).
   *
   * noInverseKalman[Z_][{x_, P_}, {A_, z_}] :=
   *   Module[{PAT, D, DiRes, DiAP, KRes, KAP},
   *    PAT = P.Transpose[A];               (* n x b *)
   *    D = Z + A.PAT;                      (* b x b *)
   *    DiRes = LinearSolve[D, z - A.x];    (* b x 1 *)
   *    KRes = PAT.DiRes;                   (* n x 1 *)
   *    DiAP = LinearSolve[D, A.P];         (* b x n *)
   *    KAP = PAT.DiAP;                     (* n x n *)
   *    {x + KRes, P - KAP}];
   */

  /* Using dgemm for P.A^T because dsymm doesn't offer a way to transpose the
     right-hand multiplicand. */

  /*
   *         A          x             z          Res
   *         n          1             1           1
   *  b / * * * * \   / * \  -->  b / * \ ,   b / * \
   *    \ * * * * / n | * |         \ * /       \ * /
   *                  | * |
   *                  \ * /
   *
   *         P           AT             PAT  
   *         n           b               b   
   *    / * * * * \   / * * \  -->    / * * \
   *  n | * * * * | n | * * |       n | * * |
   *    | * * * * |   | * * |         | * * |
   *    \ * * * * /   \ * * /         \ * * /
   *
   *       Z
   *       b
   *  b / * * \
   *    \ * * /
   *
   *         A             P           AT           APAT
   *         n             n           b              b
   *  b / * * * * \   / * * * * \   / * * \  --> b / * * \
   *    \ * * * * / n | * * * * | n | * * |        \ * * /
   *                  | * * * * |   | * * |
   *                  \ * * * * /   \ * * /
   *
   *       D         Di = D^-1   Res          DiRes
   *       b             b        1             1
   *  b / * * \ ,   b / * * \ b / * \  -->  b / * \
   *    \ * * /       \ * * /   \ * /         \ * /
   *
   *      PAT     DiRes          KRes  
   *       b        1             1   
   *    / * * \ b / * \  -->    / * \
   *  n | * * |   \ * /       n | * |
   *    | * * |                 | * |
   *    \ * * /                 \ * /
   *
   *         A     
   *         n     
   *  b / * * * * \
   *    \ * * * * /
   */

  double PAT[n * b];
  /* dgemm: http://tinyurl.com/j24npm4 */
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans,
               n, b, n,    /* m, n, k       */
               1, P, n,    /* alpha, A, pda */
               A, n,       /*        B, pdb */
               0, PAT, n); /* beta,  C, pdc */

  printf ("P.AT\n");
  for (int r = 0; r < n; ++r)
  {   for (int c = 0; c < b; ++c)
      {   printf ("%g ", PAT[c + r * b /* ncols */]);   }   }
  printf ("\n\n");

  double D[b * b];
  /* D <- A.PAT + Z (copy Z to D first) */
  cblas_dcopy (b * b, Z, 1, D, 1);
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
               b, b, n,    /* m, n, k       */
               1, A, n,    /* alpha, A, pda */
               PAT, b,     /*        B, pdb */
               1, D, b);   /* beta,  C, pdc */
  
  printf ("D\n");
  for (int r = 0; r < b; ++r)
  {   for (int c = 0; c < b; ++c)
      {   printf ("%g ", D[c + r * b]);   }   }
  printf ("\n\n");

  double Res[b * 1];
  /* Res <- (-A.x) + z (copy z to Res first)  */
  cblas_dcopy (b * 1, z, 1, Res, 1);
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
               b, 1, n,    /* m, n, k       */
               -1, A, n,   /* alpha, A, pda */
               x, 1,       /*        B, pdb */
               1, Res, 1); /* beta,  C, pdc */ 

  printf ("Res\n");
  for (int r = 0; r < b; ++r)
  {   for (int c = 0; c < 1; ++c)
      {   printf ("%g ", Res[c + r * 1]);   }   }
  printf ("\n\n");

  double DiRes[b * 1];
  double DCholesky[b * b];
  /*    DiRes = LinearSolve[D, z - A.x];    (* b x 1 *) */
  /*    copy Res to DiRes, first. */
  /*    copy D to DCholesky first. */
  /* dposv: http://goo.gl/O7gUH8 */
  cblas_dcopy (b * 1, Res, 1, DiRes,     1);
  cblas_dcopy (b * b, D,   1, DCholesky, 1);
  int result = LAPACKE_dposv (
  , 'U',
  b, /* n rows of D */
  1, /* n columns of Res */
  D,
  b,
  DiRes,
  b);

  printf ("DPOSV result %d\n\n", result);

  printf ("DiRes\n");
  for (int r = 0; r < b; ++r)
  {   for (int c = 0; c < 1; ++c)
      {   printf ("%g ", DiRes[c + r * 1]);   }   }
  printf ("\n\n");

  printf ("DCholesky\n");
  for (int r = 0; r < b; ++r)
  {   for (int c = 0; c < b; ++c)
      {   printf ("%g ", DCholesky[c + r * b]);   }   }
  printf ("\n\n");
  
  double KRes[n * 1];
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
               n, 1, b,     /* m, n, k       */
               1, PAT, 1,   /* alpha, A, pda */
               DiRes, 1,    /*        B, pdb */
               0, KRes, 1); /* beta,  C, pdc */ 

  printf ("KRes\n");
  for (int r = 0; r < n; ++r)
  {   for (int c = 0; c < 1; ++c)
      {   printf ("%g ", KRes[c + r * 1]);   }   }
  printf ("\n\n");
  
  /*  AP              = (b x n).(n x n)  -->  (b x n)
   *  DiAP            = (b x b).(b x n)  -->  (b x n)
   *  KAP  = PAT.DiAP = (n x b).(b x n)  -->  (n x n)
   */

  double AP[b * n];
  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
               b, n, n,     /* m, n, k       */
               1, A, n,     /* alpha, A, pda */
               P, n,        /*        B, pdb */
               0, AP, n);   /* beta,  C, pdc */ 

  printf ("AP\n");
  for (int r = 0; r < b; ++r)
  {   for (int c = 0; c < n; ++c)
      {   printf ("%g ", AP[c + r * n]);   }   }
  printf ("\n\n");
  
  double DiAP[b * n];
  /* DiAP = LinearSolve[D, AP];    (* b x n *) */
  /* copy AP to DiAP, first. */
  /* copy D to DCholesky first. */
  /* dposv: http://goo.gl/O7gUH8 */
  cblas_dcopy (b * n, AP, 1, DiAP,      1);
  cblas_dcopy (b * b, D,  1, DCholesky, 1);
  int result = LAPACKE_dposv (LAPACK_ROW_MAJOR, 'U', b, b, D, b, DiRes, b);

  double KAP[n * n];
  return 0;   }
