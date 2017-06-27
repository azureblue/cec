#ifndef CROSS_ENTROPY_CONTEXT_H
#define	CROSS_ENTROPY_CONTEXT_H

#include "errors.h"
#include "alloc.h"

struct cross_entropy_context
{
    res_code last_error;
    memptr_t specific_context;
};

typedef struct cross_entropy_context cross_entropy_ctx;

#endif	/* CROSS_ENTROPY_CONTEXT_H */
