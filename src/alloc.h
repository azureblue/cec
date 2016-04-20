#ifndef ALLOC_H
#define	ALLOC_H
#include <stdlib.h>

/*
 * Abstraction over memory allocation.
 */

void * m_alloc(size_t size);
void m_free(void * ptr);
void m_free_ptrs(void ** ptrs, int n);

#endif	/* ALLOC_H */
