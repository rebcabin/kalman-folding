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
#include <assert.h>

#include <Block.h>

#include <cblas.h>
#include <lapacke.h>

/*     _      _                _                   */
/*  __| |__ _| |_ _  _ _ __   | |_ _  _ _ __  ___  */
/* / _` / _` |  _| || | '  \  |  _| || | '_ \/ -_) */
/* \__,_\__,_|\__|\_,_|_|_|_|  \__|\_, | .__/\___| */
/*                                 |__/|_| */

/* I daydream and fantasize that we can abstract on numerical types. The reality
   is that a lot of this code will have secret knowledge that the underlying
   type is a floating-point scalar. Abstracting on numerical types is difficult
   even in very high-level programming languages. */

typedef double Datum;

/*           _     _                _       _      */
/*  _ __ _ _(_)_ _| |_   _ __  __ _| |_ _ _(_)_ __ */
/* | '_ \ '_| | ' \  _| | '  \/ _` |  _| '_| \ \ / */
/* | .__/_| |_|_||_\__|_|_|_|_\__,_|\__|_| |_/_\_\ */
/* |_|               |___|                         */

/* set this in options to main */
int g_verbose = 0;

void print_matrix (char * nym, const Datum * m, int rows, int cols)
{   if (! g_verbose) {
        return;   }
    printf ("%s: [\n", nym);
    for (int r = 0; r < rows; ++r)
    {   for (int c = 0; c < cols; ++c)
        {   printf ("%g, ", m[c + r * cols]);   }
        printf ("\n");   }
    printf ("],\n");   }

/*                    __     _    _      _    _       _        _                  */
/*  _ _  ___ _ _ ___ / _|___| |__| |__ _| |__| |___  | |____ _| |_ __  __ _ _ _   */
/* | ' \/ _ \ ' \___|  _/ _ \ / _` / _` | '_ \ / -_) | / / _` | | '  \/ _` | ' \  */
/* |_||_\___/_||_|  |_| \___/_\__,_\__,_|_.__/_\___| |_\_\__,_|_|_|_|_\__,_|_||_| */

/* This modifies x and P in-place. Our foldable kalman is a thin skin over this.
 * It also (conditionally) prints out all intermediate matrices for pedagogical
 * purposes. */

void kalman (int b,             /* # rows, cols, in Z; # rows in z */
             int n,             /* # rows, cols, in P; # rows in x */
             const Datum * IdN, /* n x n identity matrix */
             const Datum * Z,   /* b x b observation covariance */
             Datum * x,         /* n x 1, current state */
             Datum * P,         /* n x n, current covariance */
             const Datum * A,   /* b x n, current observation partials */
             const Datum * z    /* b x 1, current observation vector */
             ) {

    /* Transcribe the following sketch of Wolfram code (the intermediate
     * matrices are not necessary in Wolfram, but we need them in C).
     *
     * noInverseKalman[Z_][{x_, P_}, {A_, z_}] :=
     *
     *   Module[{PAT, D, Res, DiRes, KRes, AP, DiAP, KAP},
     *
     *    1. PAT    = P.Transpose[A]         (* n x b *)
     *    2. D      = Z + A.PAT              (* b x b *)
     *    3. Res    = z - A.x                (* b x 1 *)
     *    4. DiRes  = LinearSolve[D, Res]    (* b x 1 *)
     *    5. KRes   = PAT.DiRes              (* n x 1 *)
     *    6. AP     = A.P                    (* n x 1 *)
     *    7. DiAP   = LinearSolve[D, AP]     (* b x n *)
     *    8. KAP    = PAT.DiAP               (* n x n *)
     *    9. x     <- x + KRes
     *   10. P     <- P - KAP
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

    Datum PAT[n * b];
    /* dgemm: http://tinyurl.com/j24npm4 */
    /* C <-- alpha * A * B + beta * C */
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans,
                 n,          /* m (n),    # rows of A (P) */
                 b,          /* n (b),    # cols of B (AT) (post-transpose) */
                 n,          /* k (n),    # cols of A (P) == rows of B (AT post-tranpose) */
                 1, P, n,    /* alpha, A, # cols A (P,  pre-transpose)*/
                 A, n,       /*        B, # cols B (AT, pre-transpose)*/
                 0, PAT, b); /* beta,  C, # cols C */
    print_matrix ("P.AT", PAT, n, b);

    /*
     *       D                 A          PAT           Z
     *       b                 n           b            b
     *  b / * * \  <--  b / * * * * \ n / * * \  + b / * * \
     *    \ * * /         \ * * * * /   | * * |      \ * * /
     *                                  | * * |
     *                                  \ * * /
     *
     */

    Datum D[b * b];
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
    print_matrix ("D", D, b, b);

    /*
     *     Res                       A          x                 z
     *      1                        n          1                 1
     *  b / * \  <--  alpha * b / * * * * \ n / * \  + beta * b / * \
     *    \ * /                 \ * * * * /   | * |             \ * /
     *                                        | * |
     *                                        \ * /
     *
     */
    Datum Res[b * 1];
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
    print_matrix ("Res", Res, b, 1);

    /*
     *    DiRes        Di = D^-1   Res
     *      1              b        1
     *  b / * \  <--  b / * * \ b / * \
     *    \ * /         \ * * /   \ * /
     *
     */

    Datum DiRes[b * 1];
    Datum DCholesky[b * b];
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
    if (g_verbose) {
        printf ("DPOSV DiRes result %d\n\n", result);   }
    print_matrix ("DiRes",     DiRes,     b, 1);
    print_matrix ("DCholesky", DCholesky, b, b);

    /*
     *     KRes           PAT     DiRes
     *      1              b        1
     *  n / * \       n / * * \ b / * \
     *    | * |  <--    | * * |   \ * /
     *    | * |         | * * |
     *    \ * /         \ * * /
     *
     */

    Datum KRes[n * 1];
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
    print_matrix ("KRes", KRes, n, 1);

    /*
     *         AP                  A             P
     *         n                   n             n
     *  b / * * * * \  <--  b / * * * * \ n / * * * * \
     *    \ * * * * /         \ * * * * /   | * * * * |
     *                                      | * * * * |
     *                                      \ * * * * /
     *
     */

    Datum AP[b * n];
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
    print_matrix ("AP", AP, b, n);

    /*
     *        DiAP           Di = D^-1       A             P
     *         n                 b           n             n
     *  b / * * * * \  <--  b / * * \ b / * * * * \ n / * * * * \
     *    \ * * * * /         \ * * /   \ * * * * /   | * * * * |
     *                                                | * * * * |
     *                                                \ * * * * /
     *
     */

    Datum DiAP[b * n];
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
    if (g_verbose) {
        printf ("DPOSV DiAP result %d\n\n", result);   }
    print_matrix ("DiAP",      DiAP,      b, n);
    print_matrix ("DCholesky", DCholesky, b, b);

    /*
     *        KAP             PAT         DiAP
     *         n               b           n
     *  n / * * * * \  <--  / * * \ b / * * * * \
     *    | * * * * |     n | * * |   \ * * * * /
     *    | * * * * |       | * * |
     *    \ * * * * /       \ * * /
     *
     */

    Datum KAP[n * n];
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
    print_matrix ("KAP", KAP, n, n);

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
    print_matrix ("x", x, n, 1);

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
    print_matrix ("P", P, n, n); }

/*  ___  _                   _              */
/* |   \(_)_ __  ___ _ _  __(_)___ _ _  ___ */
/* | |) | | '  \/ -_) ' \(_-< / _ \ ' \(_-< */
/* |___/|_|_|_|_\___|_||_/__/_\___/_||_/__/ */

/* We give the dimensions here as constants because we define array dimensions
 * at compile time. This is an appropriate compromise between the flexibility
 * and generality of variable array dimensions and the desire to avoid heap
 * allocation as much as reasonable. If inclined to favor flexibility and
 * generality over the desire to avoid heap allocation, we would define array
 * dimensions as variables and allocate arrays themselves on the heap. If
 * inclined the other way, we would never bother storing an array dimension in a
 * variable. We split the difference: our data structures contain the array
 * dimensions as variables for easy access and for confluence with the
 * mathematical descriptions above, but we don't use the variables for
 * allocation. Instead, we must refer to the following constants when allocating
 * storage for arrays.
 *
 * We use heap allocation only at the top level and only when under warning
 * comments. It is intended to be easily replaceable by arena or stack
 * allocation.
 */

const int batch_count = 1;
const int state_count = 4;

/*    _                        _      _   _ */
/*   /_\  __ __ _  _ _ __ _  _| |__ _| |_(_)___ _ _ */
/*  / _ \/ _/ _| || | '  \ || | / _` |  _| / _ \ ' \ */
/* /_/ \_\__\__|\_,_|_|_|_\_,_|_\__,_|\__|_\___/_||_| */

typedef struct
{   int   b;
    int   n;
    Datum x[state_count];
    Datum P[state_count * state_count];   }
Accumulation, * pAccumulation;

Accumulation zeroAccumulation (void)
{   Accumulation r;
    memset ((void *)&r, 0, sizeof (Accumulation));
    return r;   }

Accumulation createAccumulation (int b_, int n_, Datum * x_, Datum * P_) {
    Accumulation r = zeroAccumulation ();
    r.b = b_;
    r.n = n_;
    assert (n_ == state_count);
    memcpy ((void *) &(r.x), (void *)x_, n_ * sizeof (Datum));
    memcpy ((void *) &(r.P), (void *)P_, n_ * n_ * sizeof (Datum));
    return r;   }

Accumulation copyAccumulation (pAccumulation pa) {
    Accumulation r;
    memcpy ((void *)&r, (void *)pa, sizeof (Accumulation));
    return r;   }

void printAccumulation (Accumulation a)
{   printf ("{b: %d, n: %d\n", a.b, a.n);
    if (! g_verbose) {
        g_verbose = 1;   }
    print_matrix ("x", a.x, a.n, 1);
    print_matrix ("P", a.P, a.n, a.n);
}

/*   ___  _                         _   _ */
/*  / _ \| |__ ___ ___ _ ___ ____ _| |_(_)___ _ _  ___ */
/* | (_) | '_ (_-</ -_) '_\ V / _` |  _| / _ \ ' \(_-< */
/*  \___/|_.__/__/\___|_|  \_/\__,_|\__|_\___/_||_/__/ */

typedef struct
{   int   n;
    Datum partials [1 * state_count];
    Datum z [batch_count * batch_count];   }
ObservationAndPartials, * pObservationAndPartials;

typedef struct
{   int count;
    int current;
    pObservationAndPartials observations_and_partials;   }
Observations;

/* private */
pObservationAndPartials allocObservationAndPartialsArray (int count_)
{   /* Don't use malloc & free in embedded apps. Use arena or stack memory. */
    pObservationAndPartials po =
        (pObservationAndPartials)
        malloc (count_ * sizeof (ObservationAndPartials));
    if (NULL == po)
    {   printf ("Failed to alloc %d observations_and_partials\n", count_);
        exit (-1);   }
    return po;   }

Observations createObservations (int count_, Datum * partials_, Datum * zs_)
{   pObservationAndPartials po = allocObservationAndPartialsArray (count_);
    for (int i = 0; i < count_; ++i) {
        po[i].n = state_count;
        memcpy ((void *) & (po[i].partials),
                & (partials_[i * state_count]),
                state_count * sizeof (Datum));
        memcpy ((void *) & (po[i].z),
                & (zs_[i]),
                sizeof (Datum));
    }
    Observations result;
    result.count   = count_;
    result.current = 0;
    result.observations_and_partials = po;
    return result;   }

void freeObservations (Observations o)
{   /* Don't use malloc & free in embedded apps. Use arena or stack memory. */
    free ((void *)o.observations_and_partials);   }

/*                       _                        */
/*  _ __ ___ ___ _  _ __| |___   ___              */
/* | '_ (_-</ -_) || / _` / _ \ |___|             */
/* | .__/__/\___|\_,_\__,_\___/                   */
/* |_|                                            */
/*              _                            _    */
/*  ___ _ ___ _(_)_ _ ___ _ _  _ __  ___ _ _| |_  */
/* / -_) ' \ V / | '_/ _ \ ' \| '  \/ -_) ' \  _| */
/* \___|_||_\_/|_|_| \___/_||_|_|_|_\___|_||_\__| */

/* In the land of real closures, free variables in the bodies of functions would
 * be "closed over," that is, copied into an environment structure, a pointer to
 * which is secretly passed as the first argument to the function (Sound
 * familiar? It's the same concept as in object-oriented programming, where a
 * pointer to the object is secretly passed to every method. In the case of
 * closures, the "object" is an environment structure created automatically by
 * the compiler by enumerating the free variables in a function body. The free
 * variables are any variables that are /not/ parameters to the function.) In
 * our case, the foldable kalman refers to two constant matrices. We'll just
 * make them static constants outside the function because their scope includes
 * the function body.
 */

static const Datum IdN[state_count * state_count] =
    {   1., 0., 0., 0.,
        0., 1., 0., 0.,
        0., 0., 1., 0.,
        0., 0., 0., 1.  };

static const Datum Z [batch_count * batch_count] = {1.};

/*    _                        _      _ */
/*   /_\  __ __ _  _ _ __ _  _| |__ _| |_ ___ _ _ */
/*  / _ \/ _/ _| || | '  \ || | / _` |  _/ _ \ '_| */
/* /_/ \_\__\__|\_,_|_|_|_\_,_|_\__,_|\__\___/_| */

typedef Accumulation (^Accumulator) (Accumulation a, ObservationAndPartials b);

Accumulator foldableKalman = ^(Accumulation a, ObservationAndPartials z) {
    /* modify a.x and a.P in-place */
    kalman (a.b, a.n, IdN, Z, a.x, a.P, z.partials, z.z);
    return a;
};

/*   __     _    _ */
/*  / _|___| |__| | */
/* |  _/ _ \ / _` | */
/* |_| \___/_\__,_| */

Accumulation fold (Accumulator f, Accumulation x0, Observations zs)
{   for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   x0 = f (x0, zs.observations_and_partials[zs.current]);   }
    return x0;   }

/*             _ */
/*  _ __  __ _(_)_ _ */
/* | '  \/ _` | | ' \ */
/* |_|_|_\__,_|_|_||_| */


int main (int argc, char ** argv)
{   Datum x[state_count * 1] =
        {0., 0., 0., 0};
    Datum P[state_count * state_count] =
        {   1000.,    0.,    0.,    0.,
               0., 1000.,    0.,    0.,
               0.,    0., 1000.,    0.,
               0.,    0.,    0., 1000.   };

    const int observations_count = 5;

    Datum partials [observations_count * state_count] =
        {    1.,  0.,  0.,  0.,
             1.,  1.,  1.,  1.,
             1., -1.,  1., -1.,
             1., -2.,  4., -8.,
             1.,  2.,  4.,  8.,   } ;
    Datum zs [observations_count * batch_count] =
        {    -2.28442,
             -4.83168,
            -10.4601,
              1.40488,
            -40.8079   };

    Accumulation initial_accumulation = createAccumulation
        (batch_count, state_count, x, P);

    Observations fu = createObservations
        (observations_count, partials, zs);

    Accumulation result = fold (foldableKalman, initial_accumulation, fu);

    printAccumulation (result);

    freeObservations (fu);

    return 0;   }
