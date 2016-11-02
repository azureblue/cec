#ifndef ALLOC_CHECK_H
#define ALLOC_CHECK_H

#include <inttypes.h>
#include <setjmp.h>
#include <stddef.h>

#define alloc_check_mem_alloc m_alloc
#define alloc_check_mem_free m_free

extern void * alloc_check_mem_alloc(size_t size);
extern void alloc_check_mem_free(void *);

struct alloc_check_context
{
    jmp_buf jmpbuf;
    const int size;
    volatile int idx;
    volatile intptr_t * ptrs;
};

/*
 * This macro implicitly declares an array of length *ptr_buf_size* (variable-length array), 
 * alloc_check_context struct and sets jmp_buf in the local scope that will be used to track 
 * and free pointers that are passed to *check_ptr* or *check_alloc*.
 */
#define checked_allocation(ptr_buf_size)                                                                                                \
    volatile intptr_t __alloc_check_ptrs[ptr_buf_size];                                                                                 \
    struct alloc_check_context __alloc_check_context = {.idx = 0, .size = ptr_buf_size, .ptrs = __alloc_check_ptrs};                    \
    if (setjmp(__alloc_check_context.jmpbuf)) return (alloc_check_free_ptrs(&__alloc_check_context), NULL);

/*
 * Macro that *alter control flow* by causing the function, inside which it has been used, 
 * to *return* NULL if *ptr* expression evaluates to NULL. Before it returns, all pointers 
 * previously passed to *check_ptr* or *check_alloc* in the context of the local scope are freed.
 *
 * The use of this macro must be preceded by *checked_allocation* in the context of the surrounding function.
 */
#define check_ptr(ptr) (alloc_check_idx_within_range(&__alloc_check_context), alloc_check_ptr(&__alloc_check_context, ptr))

/*
 * Passes alloc_check_mem_alloc(size_t_s) to check_ptr(ptr).
 * 
 * @see check_ptr
 */
#define check_alloc(size_t_s) check_ptr(alloc_check_mem_alloc(size_t_s))

#define free_checked_ptrs() alloc_check_free_ptrs(&__alloc_check_context)

void alloc_check_longjmp_clean(struct alloc_check_context *);

void alloc_check_idx_within_range(struct alloc_check_context *);

void * alloc_check_ptr(struct alloc_check_context *, void * ptr);

void alloc_check_free_ptrs(struct alloc_check_context *);

#endif /* ALLOC_CHECK_H */
