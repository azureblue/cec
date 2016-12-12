#ifndef ALLOC_H
#define	ALLOC_H
#include <stddef.h>

/*
 * Abstraction over memory allocation.
 */

typedef void* memptr_t;

memptr_t m_alloc(size_t size);
void m_free(memptr_t ptr);
void m_free_ptrs(memptr_t * ptrs, int n);

#endif	/* ALLOC_H */
