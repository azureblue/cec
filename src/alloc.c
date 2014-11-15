#include "alloc.h"

#ifdef R_ALLOC
#include <Rinternals.h>
#endif

void * m_alloc(size_t size)
{
#ifdef R_ALLOC
    return R_alloc(size / sizeof (char), 1);
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
