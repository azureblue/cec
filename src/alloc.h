#ifndef ALLOC_H
#define	ALLOC_H
#include <stddef.h>

typedef void* memptr_t;

memptr_t m_alloc(size_t size);

#define alloc(type) m_alloc(sizeof (type))
#define alloc_n(type, n) m_alloc(sizeof (type) * (n))
#define alloc_fam(struct_type, fam_type, fam_length) m_alloc(sizeof (struct_type) + sizeof (fam_type) * (fam_length))

void m_free(memptr_t ptr);
void m_free_ptrs(memptr_t * ptrs, int n);

#endif	/* ALLOC_H */
