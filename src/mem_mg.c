#include <inttypes.h>
#include "mem_mg.h"

#ifdef R_ALLOC
#include <R.h>

static memptr_t mem_alloc(size_t size) {
    return Calloc(size, char);
}

static void mem_free(memptr_t ptr) {
    Free(ptr);
}
#else
#include <stdlib.h>

static memptr_t mem_alloc(size_t size) {
    return malloc(size);
}
static void mem_free(memptr_t ptr) {
    free(ptr);
}
#endif

struct mem_mg_node {
    size_t id;
    intptr_t ptr;
    struct mem_mg_node * prev;
    struct mem_mg_node * next;
};

#define mem_node struct mem_mg_node

mem_node end_node = {.id = SIZE_MAX, .next = NULL};
mem_node start_node = {.id = 0, .prev = NULL};

static struct {
    size_t current_id;
    mem_node * start;
    mem_node * end;
    mem_fail_handler fail_handler;
} ctx = {.current_id = 1, .end = NULL};

void free_mem_mg() {
    mem_reset_state(0);
    ctx.current_id = 1;
}

void init_mem_mg(mem_fail_handler fail_handler) {
    ctx.end = &end_node;
    ctx.end->prev = &start_node;
    ctx.end->prev->next = ctx.end;
    ctx.start = &start_node;
    ctx.fail_handler = fail_handler;
}

memptr_t m_alloc(size_t size) {
    mem_node * node = mem_alloc(sizeof (mem_node));
    if (!node)
        goto fail;
    node->prev = NULL;
    node->next = NULL;
    memptr_t ptr = mem_alloc(size);
    if (!ptr)
        goto fail;
    node->id = ctx.current_id++;
    node->ptr = (intptr_t) ptr;
    node->prev = ctx.end->prev;
    node->next = ctx.end;
    node->prev->next = node;
    ctx.end->prev = node;
    return ptr;
fail:
    mem_free(node);
    ctx.fail_handler();
    return NULL;
}

mem_node * m_last_alloc_node() {
    return ctx.end->prev;
}

mem_state_id mem_track_start() {
    return ctx.end->prev->id;
}

static void free_node(mem_node *node) {
    if (!node)
        return;
    mem_free((memptr_t) node->ptr);
    if (node->prev)
        node->prev->next = node->next;
    if (node->next)
        node->next->prev = node->prev;
    mem_free(node);
}

void mem_reset_state(mem_state_id state_id) {
    mem_free_range((mem_state_range) {state_id, ctx.end->prev});
}

void mem_free_range(mem_state_range range) {
    mem_node *node = range.to_node;
    while (node->id > range.from_id) {
        mem_node * prev = node->prev;
        free_node(node);
        node = prev;
    }
}

mem_state_range mem_empty_range() {
    return (mem_state_range) {ctx.start->id, ctx.start};
}

mem_state_range mem_track_end(mem_state_id start) {
    if (m_last_alloc_node()->id <= start)
        return mem_empty_range();
    return (mem_state_range) {start, m_last_alloc_node()};
}
