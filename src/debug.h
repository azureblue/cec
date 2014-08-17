#ifndef DEBUG_H
#define	DEBUG_H

#include <stdlib.h>
#include <stdio.h>

#define malloc(T) _malloc((T))
#define free(T) _free((T))
void * _malloc(size_t size);
void _free(void *);
void reset_mallocs();
void reset_frees();
int mallocs();
int frees();
int mallocs_equals_frees();

#endif	/* DEBUG_H */

