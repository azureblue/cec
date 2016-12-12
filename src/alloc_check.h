#ifndef ALLOC_CHECK_H
#define ALLOC_CHECK_H

#include <inttypes.h>
#include <setjmp.h>
#include <stddef.h>

typedef void* memptr_t;

#define alloc_check_mem_alloc m_alloc
#define alloc_check_mem_free m_free

extern void * alloc_check_mem_alloc(size_t size);
extern void alloc_check_mem_free(void *);

typedef void (*destructor_function_t)(memptr_t);

struct alloc_check_context
{
    const int size;
    volatile int idx;
    volatile intptr_t * ptrs;
    volatile destructor_function_t * dstrs;
    jmp_buf jmpbuf;
};

#define checked_alloc(buf_size)                                                                             \
    volatile intptr_t _alloc_check_ptrs[buf_size];                                                          \
    volatile destructor_function_t _alloc_check_dstrs[buf_size];                                            \
    struct alloc_check_context _alloc_check_context = {buf_size, 0, _alloc_check_ptrs, _alloc_check_dstrs}; \
    if (setjmp(_alloc_check_context.jmpbuf)) return (alloc_check_free_ptrs(&_alloc_check_context), NULL)

#define check_ptr(ptr) (alloc_check_assert_range(&_alloc_check_context), alloc_check_ptr(&_alloc_check_context, ptr, NULL))
#define check_ptr_dstr(ptr, destr) (alloc_check_assert_range(&_alloc_check_context), alloc_check_ptr(&_alloc_check_context, ptr, destr))
#define check_alloc(struct_t) check_ptr(alloc_check_mem_alloc(sizeof (struct_t)))
#define check_alloc_n(type, n) check_ptr(alloc_check_mem_alloc(sizeof (type) * n))
#define check_alloc_fam(struct_t, fam_type, fam_length) check_ptr(alloc_check_mem_alloc(sizeof (struct_t) + sizeof (fam_type) * fam_length))

#define free_checked_ptrs() alloc_check_free_ptrs(&_alloc_check_context)

void alloc_check_longjmp_clean(struct alloc_check_context *);

void alloc_check_assert_range(struct alloc_check_context *);

memptr_t alloc_check_ptr(struct alloc_check_context *, memptr_t ptr, destructor_function_t destr);

void alloc_check_free_ptrs(struct alloc_check_context *);

#endif /* ALLOC_CHECK_H */
