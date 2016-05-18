#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Block.h>

/* Using some abbreviations lest the code be too long and harder to read. */

typedef double T, * pT;

const size_t Acmln_size = 3;
typedef T Obn, * pObn;                         /* Observation */
typedef struct s_Acmln
{   T elts[Acmln_size];   } Acmln, * pAcmln; /* Accumulation */

Acmln zeroAcmln (void)
{   Acmln r;
    for (size_t i = 0; i < Acmln_size; ++i)
    {   r.elts[i] = ((T)0);  }
    return r;   }

typedef struct s_BoundedArray_Acmlns
{   int count;
    int max;
    pAcmln * acmlns ;   } Acmlns;

void printAcmln (Acmln a)
{   printf ("{");
    for (size_t i = 0; i < Acmln_size; ++i)
    {   printf ("%lf", a.elts[i]);
        if (i < Acmln_size - 1)
        {   printf (", ");   }   }
    printf ("}\n");   }

typedef struct s_BoundedArray_Obns
{   int count;
    int current;
    Obn * obns;   } Obns;

/*private*/pObn allocObnArray (int count_)
{   /* Don't use malloc & free in embdded apps. Use arena or stack memory. */
    pObn po = (pObn) malloc (count_ * sizeof (Obn));
    if (NULL == po)
    {   printf ("Failed to alloc %d obns\n", count_);
        exit (-1);   }
    return po;   }

/* Public interface for creating arrays of observations. Part of the test
   infrastructure. */
Obns createObns (int count_, pObn pObns)
{   pObn po = allocObnArray (count_);
    memcpy ((void *)po, (void *)pObns, sizeof (Obn) * count_);
    Obns result;
    result.count   = count_;
    result.current = 0;
    result.obns    = po;
    return result;   }

/* Public interface for destroying arrays of observations. Part of the test
   infrastructure. */
void freeObns (Obns o)
{   /* Don't use malloc & free in embdded apps. Use arena or stack memory. */
    free ((void *)o.obns);   }

typedef Acmln (^Acmlr) (Acmln a, Obn b); /* Accumulator */

Acmln fold (Acmlr f, Acmln a0, Obns zs)
{   for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   a0 = f (a0, zs.obns[zs.current]);   }
    return a0;   }

pAcmln foldList (Acmlr f, Acmln a0, Obns zs)
{   return NULL;   }

int main (int argc, char ** argv)
{   Obn tmp[3] = {55, 89, 144};
    Obns zs = createObns(3, tmp);
    Acmln x0 = zeroAcmln ();
    Acmln result = fold (^(Acmln a, Obn z)
                         {   T var = a.elts[0];
                             T x   = a.elts[1];
                             T n   = a.elts[2];

                             T K = 1.0 / (1.0 + n);
                             T x2 = x + K * (z - x);
                             T ssr2 = (n - 1.0) * var + K * n * (z - x) * (z - x);

                             Acmln r;
                             r.elts[0] = ssr2 / (n > 1.0 ? n : 1.0);
                             r.elts[1] = x2;
                             r.elts[2] = n + 1.0;
                             return r;},
                         x0,
                         zs);
    printAcmln (result);
    freeObns (zs);
    return 0;   }
