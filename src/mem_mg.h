#ifndef CEC_CMEM_MG_H
#define CEC_CMEM_MG_H
#include <inttypes.h>
#include <setjmp.h>
#include <stddef.h>

typedef void* memptr_t;
typedef void (*mem_fail_handler)();
typedef size_t m_state;

enum mem_mg_init_res {
    OK = 0,
    ALREADY_INITIALIZED = 1,
    FAILED = -1
};

enum mem_mg_init_res init_mem_mg(mem_fail_handler fail_handler);
void free_mem_mg();
memptr_t m_alloc(size_t size);
m_state m_current_state();
void m_reset_state(m_state);

#endif //CEC_CMEM_MG_H
