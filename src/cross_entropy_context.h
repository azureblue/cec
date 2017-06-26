#ifndef CROSS_ENTROPY_CONTEXT_H
#define	CROSS_ENTROPY_CONTEXT_H

#include "errors.h"
#include "alloc.h"

enum density_family
{
    GIVEN_COVARIANCE = 0, FIXED_R = 1, SPHERICAL = 2, DIAGONAL = 3, FIXEDEIGENVALUES = 4, ALL = 5
};

struct cross_entropy_context
{
    res_code last_error;
    memptr_t specific_context;
};

typedef struct cross_entropy_context cross_entropy_ctx;

#endif	/* CROSS_ENTROPY_CONTEXT_H */
