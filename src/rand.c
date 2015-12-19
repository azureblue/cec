#include <stdlib.h>
#include "rand.h"
#ifdef R_RAND
#include<R_ext/Random.h>
#endif

void cec_rand_init()
{
#ifdef R_RAND
    GetRNGstate();
#endif
}

void cec_rand_end()
{
#ifdef R_RAND
    PutRNGstate();
#endif
}

double cec_rand()
{
#ifdef R_RAND
    return unif_rand();
#else
    return (double) rand() / (double) RAND_MAX;
#endif
}
