#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Block.h>
typedef double T;

const size_t Accumulation_size = 2;
typedef struct s_Accumulation
{   T elements[Accumulation_size];   } Accumulation;

Accumulation zeroAccumulation (void)
{   Accumulation r;
    memset ((void *)r.elements, 0, Accumulation_size * sizeof (T));
    return r;   }

void printAccumulation (Accumulation a)
{    printf ("{%lf, %lf}\n", a.elements[0], a.elements[1]);   }
typedef T Observation, * pObservation;
typedef struct s_BoundedArray_Observations
{   int count;
    int current;
    pObservation observations;   } Observations;

/*private*/pObservation allocObservationArray (int count_)
{   /* Don't use malloc & free in embedded apps. Use arena or stack memory. */
    pObservation po = (pObservation) malloc (count_ * sizeof (Observation));
    if (NULL == po)
    {   printf ("Failed to alloc %d observations\n", count_);
        exit (-1);   }
    return po;   }

Observations createObservations (int count_, pObservation pObservations)
{   pObservation po = allocObservationArray (count_);
    memcpy ((void *)po, (void *)pObservations, sizeof (Observation) * count_);
    Observations result;
    result.count   = count_;
    result.current = 0;
    result.observations    = po;
    return result;   }

void freeObservations (Observations o)
{   /* Don't use malloc & free in embedded apps. Use arena or stack memory. */
    free ((void *)o.observations);   }
typedef Accumulation (^Accumulator) (Accumulation a, Observation b);
Accumulation fold (Accumulator f, Accumulation x0, Observations zs)
{   for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   x0 = f (x0, zs.observations[zs.current]);   }
    return x0;   }
int main (int argc, char ** argv)
{   Accumulator cume = ^(Accumulation a, Observation z)
        {   /* unpack inputs */
            T x = a.elements[0];
            T n = a.elements[1];
            /* compute gain */
            T K = 1.0 / (1.0 + n);
            /* busines logic, and packing results */
            Accumulation r;
            r.elements[0] = x + K * (z - x);
            r.elements[1] = 1.0 + n;

            return r;   };

    Accumulation x0 = zeroAccumulation ();

    Observation tmp[3] = {55, 89, 144};
    Observations zs = createObservations(3, tmp);

    Accumulation result = fold (cume, x0, zs);

    printAccumulation (result);
    return 0;   }
