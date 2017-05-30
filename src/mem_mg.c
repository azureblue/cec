#include <stdbool.h>
#include "mem_mg.h"

#define DEFAULT_START_SIZE 10

#ifdef R_ALLOC
#include <R.h>
static memptr_t mem_alloc(size_t size) {
    return Calloc(size, char);
}

static memptr_t mem_realloc(memptr_t ptr, size_t size) {
    return Realloc(ptr, size, char);
}
static void mem_free(memptr_t ptr) {
    Free(ptr);
}

#else
#include <stdlib.h>

static memptr_t mem_alloc(size_t size) {
    return malloc(size);
}

static memptr_t mem_realloc(memptr_t ptr, size_t size) {
    return realloc(ptr, size);
}

static void mem_free(memptr_t ptr) {
    free(ptr);
}
#endif

static struct {
    size_t size;
    size_t idx;
    intptr_t **ptrs;
    mem_fail_handler fail_handler;
} mem_mg_ctx = {.size = 0, .idx = 0, .ptrs = NULL};

static void handle_fail(mem_fail_handler fail_handler) {
    fail_handler();
}

void free_mem_mg() {
    if (mem_mg_ctx.size == 0)
        return;
    m_reset_state(0);
    mem_free(mem_mg_ctx.ptrs);
    mem_mg_ctx.size = 0;
    mem_mg_ctx.idx = 0;
    mem_mg_ctx.ptrs = NULL;
}

enum mem_mg_init_res init_mem_mg(mem_fail_handler fail_handler) {
    if (mem_mg_ctx.size > 0)
        return ALREADY_INITIALIZED;
    intptr_t **ptrs = mem_alloc(DEFAULT_START_SIZE * sizeof(intptr_t));
    if (!ptrs) {
        handle_fail(fail_handler);
        return FAILED;
    }
    mem_mg_ctx.fail_handler = fail_handler;
    mem_mg_ctx.ptrs = ptrs;
    mem_mg_ctx.size = DEFAULT_START_SIZE;
    mem_mg_ctx.idx = 0;
    return OK;
}

static bool mem_mg_grow() {
    size_t nsize = mem_mg_ctx.size * 2;
    intptr_t **ptrs = mem_realloc((memptr_t) mem_mg_ctx.ptrs, sizeof(intptr_t) * nsize);
    if (!ptrs)
        return false;
    mem_mg_ctx.ptrs = ptrs;
    mem_mg_ctx.size = nsize;
    return true;
}

memptr_t m_alloc(size_t size) {
    if (mem_mg_ctx.idx == mem_mg_ctx.size && !mem_mg_grow())
        goto fail;
    memptr_t ptr = mem_alloc(size);
    if (!ptr)
        goto fail;
    mem_mg_ctx.ptrs[mem_mg_ctx.idx++] = ptr;
    return ptr;
fail:
    mem_mg_ctx.fail_handler();
    return NULL;
}

m_state m_current_state() {
    return mem_mg_ctx.idx;
}

void m_reset_state(m_state idx) {
    while (mem_mg_ctx.idx-- > idx)
        mem_free(mem_mg_ctx.ptrs[mem_mg_ctx.idx]);
    mem_mg_ctx.idx = idx;
}

void m_clear_states(m_state a, m_state b) {
    m_state tmp;
    if (a > b) {
        tmp = a;
        a = b;
        b = tmp;
    }

    for (m_state i = a; i < b; i++) {
        mem_free(mem_mg_ctx.ptrs[i]);
        mem_mg_ctx.ptrs[i] = NULL;
    }        
}
