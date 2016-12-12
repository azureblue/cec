#include "alloc_check.h"

void alloc_check_longjmp_clean(struct alloc_check_context * ctx)
{
    longjmp(ctx->jmpbuf, 1);
}

void alloc_check_assert_range(struct alloc_check_context * ctx)
{
    if (ctx->idx >= ctx->size)
        alloc_check_longjmp_clean(ctx);
}

memptr_t alloc_check_ptr(struct alloc_check_context * ctx, memptr_t ptr, destructor_function_t dstr)
{
    if (!ptr)
        alloc_check_longjmp_clean(ctx);
    
    ctx->ptrs[ctx->idx] = (intptr_t) ptr;
    ctx->dstrs[ctx->idx] = dstr;
    
    ctx->idx++;
    
    return ptr;
}

void alloc_check_free_ptrs(struct alloc_check_context * ctx)
{
    for (int i = 0; i < ctx->idx; i++)
    {
        destructor_function_t dstr = ctx->dstrs[i];
        memptr_t ptr = (memptr_t) ctx->ptrs[i];
        
        if (dstr != NULL)
            dstr(ptr);
        else
            alloc_check_mem_free(ptr);
    }
    ctx->idx = 0;
}
