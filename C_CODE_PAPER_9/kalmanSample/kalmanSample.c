/*
  Copyright 2016 Brian C. Beckman

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/
/* This is an educational example only, not suitable for real applications.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>

void printm (char * nym, double * m, int rows, int cols)
{   printf ("%s\n", nym);
    for (int r = 0; r < rows; ++r)
    {   for (int c = 0; c < cols; ++c)
        {   printf ("%g ", m[c + r * cols]);   }
        printf ("\n");   }
    printf ("\n");   }

void kalman (int b,        /* # rows, cols, in Z; # rows in z */
             int n,        /* # rows, cols, in P; # rows in x */
             double * IdN, /* n x n identity matrix */
             double * Z,   /* b x b observation covariance */
             double * x,   /* n x 1, current state */
             double * P,   /* n x n, current covariance */
             double * A,   /* b x n, current observation partials */
             double * z    /* b x 1, current observation vector */
             ) {

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


    /* Use dgemm for P.A^T because dsymm doesn't offer a way to transpose the
       right-hand multiplicand. */

    /*
     *      PAT                P           AT
     *       b                 n           b
     *    / * * \         / * * * * \   / * * \
     *  n | * * |  <--  n | * * * * | n | * * |
     *    | * * |         | * * * * |   | * * |
     *    \ * * /         \ * * * * /   \ * * /
     *
     */

    double PAT[n * b];
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans,
                 n,          /* m (n),    # rows of A (P) */
                 b,          /* n (b),    # cols of B (AT) (post-transpose) */
                 n,          /* k (n),    # cols of A (P) == rows of B (AT post-tranpose) */
                 1, P, n,    /* alpha, A, # cols A (P,  pre-transpose)*/
                 A, n,       /*        B, # cols B (AT, pre-transpose)*/
                 0, PAT, b); /* beta,  C, # cols C */
    printm ("P.AT", PAT, n, b);

    /*
     *       D                 A          PAT           Z
     *       b                 n           b            b
     *  b / * * \  <--  b / * * * * \ n / * * \  + b / * * \
     *    \ * * /         \ * * * * /   | * * |      \ * * /
     *                                  | * * |
     *                                  \ * * /
     *
     */

    double D[b * b];
    /* D <- A.PAT + Z (copy Z to D first) */
    cblas_dcopy (b * b, Z, 1, D, 1);
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 b,          /* m (b),          # rows of A (A) */
                 b,          /* n (b),          # cols of B (PAT) */
                 n,          /* k (n),          # cols of A (A) == rows of B (PAT) */
                 1, A, n,    /* alpha, A (A),   # cols A (A) */
                 PAT, b,     /*        B (PAT), # cols B (PAT)*/
                 1, D, b);   /* beta,  C (Z),   # cols C (D) */
    printm ("D", D, b, b);

    /*
     *     Res                       A          x                 z
     *      1                        n          1                 1
     *  b / * \  <--  alpha * b / * * * * \ n / * \  + beta * b / * \
     *    \ * /                 \ * * * * /   | * |             \ * /
     *                                        | * |
     *                                        \ * /
     *
     */
    double Res[b * 1];
    /* Res <- (-A.x) + z (copy z to Res first)  */
    cblas_dcopy (b * 1, z, 1, Res, 1);
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 b,          /* m (b),        # rows of A (A) */
                 1,          /* n (1),        # cols of B (x) */
                 n,          /* k (n),        # cols of A (A) == rows of B (x) */
                 -1, A, n,   /* alpha, A (A), # cols A (A) */
                 x, 1,       /*        B (x), # cols B (x) */
                 1, Res, 1); /* beta,  C (z), # cols C (Res) */
    printm ("Res", Res, b, 1);

    /*
     *    DiRes        Di = D^-1   Res
     *      1              b        1
     *  b / * \  <--  b / * * \ b / * \
     *    \ * /         \ * * /   \ * /
     *
     */

    double DiRes[b * 1];
    double DCholesky[b * b];
    /* DiRes = LinearSolve[D, z - A.x];    (* b x 1 *) */
    /* copy Res to DiRes, first. */
    /* copy D to DCholesky first. */
    /* dposv: http://goo.gl/O7gUH8 */
    cblas_dcopy (b * 1, Res, 1, DiRes,     1);
    cblas_dcopy (b * b, D,   1, DCholesky, 1);
    int result = LAPACKE_dposv (LAPACK_ROW_MAJOR, 'U',
                                b,         /* NEQS: # rows of D */
                                1,         /* NRHS: # columns of z - A.x == Res */
                                DCholesky, /* DCholesky starts as D */
                                b,         /* PDA D */
                                DiRes,     /* output buffer */
                                b);        /* PDA DiRes */
    printf ("DPOSV DiRes result %d\n\n", result);
    printm ("DiRes",     DiRes,     b, 1);
    printm ("DCholesky", DCholesky, b, b);

    /*
     *     KRes           PAT     DiRes
     *      1              b        1
     *  n / * \       n / * * \ b / * \
     *    | * |  <--    | * * |   \ * /
     *    | * |         | * * |
     *    \ * /         \ * * /
     *
     */

    double KRes[n * 1];
    /* KRes <-- PAT.DiRes */
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 n,           /* m (n),            # rows of A (PAT) */
                 1,           /* n (1),            # cols of B (DiRes) */
                 b,           /* k (b),            # cols of A (PAT) == # rows of B (DiRes) */
                 1, PAT, b,   /* alpha, A (PAT),   # cols A (PAT) */
                 DiRes, 1,    /*        B (DiRes), # cols B (DiRes) */
                 0, KRes, 1); /* beta,  C (KRes),  # cols C (KRes) */
    printm ("KRes", KRes, n, 1);

    /*
     *         AP                  A             P
     *         n                   n             n
     *  b / * * * * \  <--  b / * * * * \ n / * * * * \
     *    \ * * * * /         \ * * * * /   | * * * * |
     *                                      | * * * * |
     *                                      \ * * * * /
     *
     */

    double AP[b * n];
    /* AP <-- A.P */
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 b,           /* m (b),          # rows of A (A) */
                 n,           /* n (n),          # cols of B (P) */
                 n,           /* k (n),          # cols of A (A) == # rows of B (P) */
                 1, A, n,     /* alpha, A (A),   # cols A (PAT) */
                 P, n,        /*        B (P),   # cols B (DiRes) */
                 0, AP, n);   /* beta,  C (AP),  # cols C (KRes) */
    printm ("AP", AP, b, n);

    /*
     *        DiAP           Di = D^-1       A             P
     *         n                 b           n             n
     *  b / * * * * \  <--  b / * * \ b / * * * * \ n / * * * * \
     *    \ * * * * /         \ * * /   \ * * * * /   | * * * * |
     *                                                | * * * * |
     *                                                \ * * * * /
     *
     */

    double DiAP[b * n];
    /* DiAP = LinearSolve[D, AP];    (* b x n *) */
    /* copy AP to DiAP, first. */
    /* copy D to DCholesky first. */
    /* dposv: http://goo.gl/O7gUH8 */
    cblas_dcopy (b * n, AP, 1, DiAP,      1);
    cblas_dcopy (b * b, D,  1, DCholesky, 1);
    result = LAPACKE_dposv (LAPACK_ROW_MAJOR, 'U',
                            b,         /* NEQS: # rows of D */
                            n,         /* NRHS: # columns of z - A.x == Res */
                            DCholesky, /* DCholesky starts as D */
                            b,         /* PDA D */
                            DiAP,      /* output buffer */
                            n);        /* PDA DiRes */
    printf ("DPOSV DiAP result %d\n\n", result);
    printm ("DiAP",      DiAP,      b, n);
    printm ("DCholesky", DCholesky, b, b);

    /*
     *        KAP             PAT         DiAP
     *         n               b           n
     *  n / * * * * \  <--  / * * \ b / * * * * \
     *    | * * * * |     n | * * |   \ * * * * /
     *    | * * * * |       | * * |
     *    \ * * * * /       \ * * /
     *
     */

    double KAP[n * n];
    /* KAP <-- PAT.DiAP */
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 n,             /* m (n),           # rows of A (PAT) */
                 n,             /* n (n),           # cols of B (DiAP) */
                 b,             /* k (b),           # cols of A (PAT) == # rows of B (DiAP) */
                 1, PAT, b,     /* alpha, A (PAT),  # cols A (PAT) */
                 DiAP, n,       /*        B (Diap), # cols B (DiRes) */
                 0, KAP, n);    /* beta,  C (KAP),  # cols C (KAP) */
    printm ("KAP", KAP, n, n);

    /*
     *      x                       Id          x                 KRes
     *      1                        n          1                  1
     *  n / * \  <--  alpha * n / * * * * \ n / * \  +  beta * n / * \
     *    | * |                 | * * * * |   | * |              | * |
     *    | * |                 | * * * * |   | * |              | * |
     *    \ * /                 \ * * * * /   \ * /              \ * /
     *
     */

    /* x <-- alpha * IdN[n] * KRes + beta * x */
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 n,           /* m (n),           # rows of A (Id) */
                 1,           /* n (1),           # cols of B (x) */
                 n,           /* k (n),           # cols of A (Id) == rows of B (x) */
                 1, IdN, n,   /* alpha, A (Id),   # cols A */
                 x, 1,        /*        B (x),    # cols B */
                 1, KRes, 1); /* beta,  C (Kres), # cols C (new x) */
    cblas_dcopy (n * 1, KRes, 1, x, 1);
    printm ("x", x, n, 1);

    /*
     *         P                          Id            KAP                       P
     *         n                           n             n                        n
     *  n / * * * * \  <--  alpha * n / * * * * \ n / * * * * \  +  beta * n / * * * * \
     *    | * * * * |                 | * * * * |   | * * * * |              | * * * * |
     *    | * * * * |                 | * * * * |   | * * * * |              | * * * * |
     *    \ * * * * /                 \ * * * * /   \ * * * * /              \ * * * * /
     *
     */

    /* P <-- P - KAP == - IdN[n] * KAP  + P */
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans,
                 n,           /* m (n),         # rows of A (Id) */
                 n,           /* n (n),         # cols of B (KAP) */
                 n,           /* k (n),         # cols of A (Id) == rows of B (KAP) */
                 -1, IdN, n,  /* alpha, A (Id), # cols A */
                 KAP, n,      /*        B (x),  # cols B */
                 1, P, n);    /* beta,  C (P),  # cols C (new P) */
    printm ("P", P, n, n); }

int main (int argc, char ** argv)
{   const int    b = 1;
    const int    n = 4;

    double IdN[n * n] = { 1., 0., 0., 0.,
                          0., 1., 0., 0.,
                          0., 0., 1., 0.,
                          0., 0., 0., 1. };


    double Z[b * b] = {1.};

    double x[n * 1] = {0., 0., 0., 0};
    double P[n * n] = {1000.,    0.,    0.,    0.,
                       0., 1000.,    0.,    0.,
                       0.,    0., 1000.,    0.,
                       0.,    0.,    0., 1000. };

    double A[b * n] = {1., 0., 0., 0};
    double z[b] = {-2.28442};

    kalman (b, n, IdN, Z, x, P, A, z);

    A[0] = 1;
    A[1] = 1;
    A[2] = 1;
    A[3] = 1;

    z[0] = -4.83168;

    kalman (b, n, IdN, Z, x, P, A, z);

    A[0] = 1;
    A[1] = -1;
    A[2] = 1;
    A[3] = -1;

    z[0] = -10.4601;

    kalman (b, n, IdN, Z, x, P, A, z);

    A[0] = 1;
    A[1] = -2;
    A[2] = 4;
    A[3] = -8;

    z[0] = 1.40488;

    kalman (b, n, IdN, Z, x, P, A, z);

    A[0] = 1;
    A[1] = 2;
    A[2] = 4;
    A[3] = 8;

    z[0] = -40.8079;

    kalman (b, n, IdN, Z, x, P, A, z);

    return 0;   }
