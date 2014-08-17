#include "debug.h"
#include <stdio.h>
#include <R.h>

#undef malloc
#undef free

static int ml = 0;
static int fs = 0;

void * _malloc(size_t size)
{
	ml++;
	return malloc(size);
}

void _free(void * ptr)
{
	if (ptr != NULL)
		fs++;
	free(ptr);
}

#define malloc(T) _malloc((T))
#define free(T) _free((T))

void reset_mallocs()
{
	ml = 0;
}

void reset_frees()
{
	fs = 0;
}

int mallocs()
{
	return ml;
}

int frees()
{
	return fs;
}

int mallocs_equals_frees()
{
	return ml == fs;
}
