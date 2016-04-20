#include "alloc.h"
#include <stdlib.h>

#ifdef R_ALLOC
#include <Rinternals.h>
#endif

void * m_alloc(size_t size)
{
#ifdef R_ALLOC
    return R_alloc(size, 1);
#else
    return malloc(size);
#endif    
}

void m_free(void * ptr)
{
#ifdef R_ALLOC
    return;
#else
    free(ptr);
#endif
}

void m_free_ptrs(void ** ptrs, int n)
{
    for (int i = 0; i < n; i++)
        m_free(ptrs[i]);
}

