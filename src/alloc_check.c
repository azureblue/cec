#include "alloc_check.h"

void alloc_check_longjmp_clean(struct alloc_check_context * context)
{
    longjmp(context->jmpbuf, 1);
}

void alloc_check_idx_within_range(struct alloc_check_context * context)
{
    if (context->idx >= context->size)
        alloc_check_longjmp_clean(context);
}

void * alloc_check_ptr(struct alloc_check_context * context, void * ptr)
{
    context->ptrs[context->idx++] = (intptr_t) ptr;
    if (!ptr)
        alloc_check_longjmp_clean(context);
    
    return ptr;
}

void alloc_check_free_ptrs(struct alloc_check_context * context)
{
    for (int i = 0; i < context->idx; i++)
        alloc_check_mem_free((void *) context->ptrs[i]);

    context->idx = 0;
}
