#ifndef ALLOC_H
#define	ALLOC_H
#include "mem_mg.h"

#define alloc(type) m_alloc(sizeof (type))
#define alloc_n(type, n) m_alloc(sizeof (type) * (n))
#define alloc_fam(struct_type, fam_type, fam_length) m_alloc(sizeof (struct_type) + sizeof (fam_type) * (fam_length))

#endif	/* ALLOC_H */
