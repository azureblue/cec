#ifndef MEM_MG_H
#define MEM_MG_H
#include <stddef.h>

typedef void* memptr_t;
typedef void (*mem_fail_handler)();
typedef size_t mem_state_id;

struct mem_mg_state_range {
    size_t from_id;
    memptr_t to_node;
};
typedef struct mem_mg_state_range mem_state_range;

memptr_t m_alloc(size_t size);

void init_mem_mg(mem_fail_handler fail_handler);
void free_mem_mg();

void mem_free_range(mem_state_range);
void mem_reset_state(mem_state_id);

mem_state_id mem_track_start();
mem_state_range mem_track_end(mem_state_id start);
mem_state_range mem_empty_range();

#define mem_track(mem_alloc_exp, m_range_ptr) {   \
    mem_state_id start = mem_track_start();       \
    mem_alloc_exp;                                \
    *(m_range_ptr) = mem_track_end(start);        \
}

#endif //MEM_MG_H
