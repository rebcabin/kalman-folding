#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Block.h>

/* Using some abbreviations lest the code be too long and harder to read. */

typedef double T;
typedef T   Obn, * pObn;  /* Observation */
typedef T   Acmln;        /* Accumulation */

typedef struct s_BoundedArray
{   int count;
    int current;
    pObn obns;   } Obns;

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

typedef Acmln (^ Acmlr) (Acmln a, Obn b); /* Accumulator */

Acmlr mul = (^(Acmln a, Obn z) {return ((Acmln)(a * z));});

Acmln fold (Acmlr f, Acmln x0, Obns zs)
{   for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   x0 = f (x0, zs.obns[zs.current]);   }
    return x0;   }

int main (int argc, char ** argv)
{   Obn tmp[3] = {55, 89, 144};
    Obns zs = createObns(3, tmp);

    printf ("%lf\n", fold (mul, 1, zs));
    printf ("%lf\n", fold ( (^(Acmln a, Obn b) {return a + b;}), 0, zs));
    printf ("%lf\n", fold ( (^(Acmln a, Obn b) {return a + 1;}), 0, zs));

    freeObns (zs);

    return 0;   }
