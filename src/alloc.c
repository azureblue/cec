#include "alloc.h"

#ifdef R_ALLOC
#include <Rinternals.h>

void * m_alloc(size_t size)
{
	return R_alloc(size / sizeof(char), 1);
    }
    
void m_free(void * ptr)
{
    return;
}

#else

void * m_alloc(size_t size)
{
    return malloc(size);
}

void m_free(void * ptr)
{
    free(ptr);
}

#endif
