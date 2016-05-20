#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Block.h>

typedef double T;

const size_t Accumulation_size = 3;
typedef T Observation, * pObservation;
typedef struct s_Accumulation
{   T elements[Accumulation_size];   } Accumulation, * pAccumulation;

Accumulation zeroAccumulation (void)
{   Accumulation r;
    memset ((void *)r.elements, 0, Accumulation_size * sizeof (T));
    return r;   }

void printAccumulation (Accumulation a)
{   printf ("{");
    for (size_t i = 0; i < Accumulation_size; ++i)
    {   printf ("%lf", a.elements[i]);
        if (i < Accumulation_size - 1)
        {   printf (", ");   }   }
    printf ("}\n");   }

/* This structure must manage memory because we do not statically know how many
   of them there are. */
typedef struct s_BoundedArray_Accumulations
{   int count;
    int max;
    pAccumulation accumulations ;   } Accumulations;

Accumulation lastAccumulations (Accumulations as)
{   if (0 == as.count)
    {   printf ("Attempt to pull non-existent element\n");
        exit (-4);   }
    return as.accumulations[as.count - 1];   }

Accumulations appendAccumulations (Accumulations as, Accumulation a)
{   Accumulations result = as;
    if (result.count + 1 > result.max)
    {   /* Double the storage. */
        int new_max = 2 * result.max;
        /* Don't use malloc & free in embdded apps. Use arena or stack memory. */
        pAccumulation new = (pAccumulation) malloc (sizeof (Accumulation) * new_max);
        if (NULL == new)
        {   printf ("Failed to alloc %d Accumulations\n", new_max);
            exit (-2);   }
        if (result.count != result.max)
        {   printf ("Internal bugcheck\n");
            exit (-3);   }
        memset ((void *)new, 0, new_max * sizeof (Accumulation));
        memcpy ((void *)new, (void *)result.accumulations, (sizeof (Accumulation) * result.max));
        free ((void *) result.accumulations);
        result.accumulations = new; 
        result.max = new_max;   }
    result.accumulations[result.count] = a;
    ++ result.count;
    return result;
}

Accumulations createAccumulations (void)
{   Accumulations result;
    const int init_size = 4;
    result.max = init_size;
    result.count = 0;
    result.accumulations = (pAccumulation) malloc (sizeof (Accumulation) * init_size);
    memset ((void *)result.accumulations, 0, sizeof (Accumulation) * init_size);
    return result;   }

void freeAccumulations (Accumulations as)
{   memset ((void *) as.accumulations, 0, (sizeof (Accumulation) * as.count)); 
    free ((void *) as.accumulations);   }

void printAccumulations (Accumulations as)
{   for (int j = 0; j < as.count; ++j )
    {   printAccumulation (as.accumulations[j]);   }   }

typedef struct s_BoundedArray_Observations
{   int count;
    int current;
    Observation * obns;   } Observations;

/*private*/pObservation allocObservationArray (int count_)
{   /* Don't use malloc & free in embdded apps. Use arena or stack memory. */
    pObservation po = (pObservation) malloc (count_ * sizeof (Observation));
    if (NULL == po)
    {   printf ("Failed to alloc %d obns\n", count_);
        exit (-1);   }
    return po;   }

/* Public interface for creating arrays of observations. Part of the test
   infrastructure. */
Observations createObservations (int count_, pObservation pObservations)
{   pObservation po = allocObservationArray (count_);
    memcpy ((void *)po, (void *)pObservations, sizeof (Observation) * count_);
    Observations result;
    result.count   = count_;
    result.current = 0;
    result.obns    = po;
    return result;   }

/* Public interface for destroying arrays of observations. Part of the test
   infrastructure. */
void freeObservations (Observations o)
{   /* Don't use malloc & free in embdded apps. Use arena or stack memory. */
    free ((void *)o.obns);   }

typedef Accumulation (^Accumulator) (Accumulation a, Observation b); /* Accumulator */

Accumulation fold (Accumulator f, Accumulation a0, Observations zs)
{   for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   a0 = f (a0, zs.obns[zs.current]);   }
    return a0;   }

Accumulations foldList (Accumulator f, Accumulation a0, Observations zs)
{   Accumulations result = createAccumulations ();
    result = appendAccumulations (result, a0);
    for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   result = appendAccumulations (result, f(lastAccumulations(result), zs.obns[zs.current]));   }
        return result;   }

int main (int argc, char ** argv)
{   Observation tmp[3] = {55, 89, 144};
    Observations zs = createObservations(3, tmp);
    Accumulation x0 = zeroAccumulation ();
    Accumulator cume = ^(Accumulation a, Observation z)
        {   T var = a.elements[0];
            T x   = a.elements[1];
            T n   = a.elements[2];

            T K = 1.0 / (1.0 + n);
            T x2 = x + K * (z - x);
            T ssr2 = (n - 1.0) * var + K * n * (z - x) * (z - x);

            Accumulation r;
            r.elements[0] = ssr2 / (n > 1.0 ? n : 1.0);
            r.elements[1] = x2;
            r.elements[2] = n + 1.0;
            return r;   };

    Accumulation result = fold (cume, x0, zs);
    printAccumulation (result);

    Accumulations results = foldList (cume, x0, zs);
    printAccumulations (results);

    freeAccumulations (results);
    freeObservations (zs);
    return 0;   }
