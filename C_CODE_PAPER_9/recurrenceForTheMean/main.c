#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Block.h>

/* Using some abbreviations lest the code be too long and harder to read. */

typedef double T, * pT;

const size_t Acmln_size = 2;
typedef T Obn, * pObn;                         /* Observation */
typedef struct s_Acmln
{   T acmlns[Acmln_size];   } Acmln, * pAcmln; /* Accumulation */
/*
pAcmln cpyAcmln (pAcmln src)
{   pAcmln dest = (pAcmln) malloc (sizeof (Acmln));
    memcpy ((void *)dest, (void *)src, sizeof (Acmln));
    return dest;   }

void freeAcmln (pAcmln pa)
{    free ((void *) pa);   }
*/
void printAcmln (Acmln a)
{    printf ("{%lf, %lf}\n", a.acmlns[0], a.acmlns[1]);   }

typedef struct s_BoundedArray
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

Obns createObns (int count_, pObn pObns)
{   pObn po = allocObnArray (count_);
    memcpy ((void *)po, (void *)pObns, sizeof (Obn) * count_);
    Obns result;
    result.count   = count_;
    result.current = 0;
    result.obns    = po;
    return result;   }

void freeObns (Obns o)
{   /* Don't use malloc & free in embdded apps. Use arena or stack memory. */
    free ((void *)o.obns);   }

typedef Acmln (^Acmlr) (Acmln a, Obn b); /* Accumulator */

// Acmlr cume = (^(Acmln a, Obn z) {return ((Acmln)(a * z));});

Acmln fold (Acmlr f, Acmln x0, Obns zs)
{   for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   x0 = f (x0, zs.obns[zs.current]);   }
    return x0;   }

int main (int argc, char ** argv)
{   Obn tmp[3] = {55, 89, 144};
    Obns zs = createObns(3, tmp);

    Acmln x0;
    x0.acmlns[0] = 0.0;
    x0.acmlns[1] = 0.0;

    Acmln result = fold (^(Acmln x, Obn z)
                         {   T K = 1.0 / (1.0 + x.acmlns[1]);
                             Acmln r;
                             r.acmlns[0] = x.acmlns[0] + K * (z - x.acmlns[0]);
                             r.acmlns[1] = 1.0 + x.acmlns[1];
                             return r;},
                         x0,
                         zs);
    printAcmln (result);
    return 0;   }
