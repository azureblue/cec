#ifndef CHECKED_ALLOCATOR_H
#define CHECKED_ALLOCATOR_H

#include "alloc.h"

#define checked_allocator_size 100

/*
 * This macro implicitly declares an array and an integer in
 * the local scope that will be used to track and free pointers
*  that are passed to the subsequent uses of the *check* macro.
 */
#define checked_allocation void * __checked_allocator_ptrs[checked_allocator_size]; int __checked_allocator_idx = 0;

/*
 * Frees pointers thet were previously passed to the *check* macro.
 */
#define free_checked_pointers() m_free_ptrs(__checked_allocator_ptrs, __checked_allocator_idx)

#define free_checked_pointers_and_return_null() {free_checked_pointers(); return NULL;}

/*
 * Macro that *alters control flow* by causing the function, inside
 * which it has been used, to *return* NULL if the *pointer* expression
 * evaluates to NULL. Before the return it frees all pointers that
 * were previously passed to the *check* in the context of the local scope.
 *
 * This macro should be used as the rightmost expression of an assignment.
 *
 * The use of this macro must be preceded by the use of *checked_allocation;* 
 * in the context of the surrounding function.
 */
#define check(pointer) \
__checked_allocator_ptrs[__checked_allocator_idx++] = pointer; \
if (__checked_allocator_idx == checked_allocator_size) free_checked_pointers_and_return_null(); \
if (!__checked_allocator_ptrs[__checked_allocator_idx - 1]) free_checked_pointers_and_return_null();

#endif /* CHECKED_ALLOCATOR_H */
