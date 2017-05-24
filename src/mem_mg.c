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
} ctx = {.size = 0, .idx = 0, .ptrs = NULL};

static void handle_fail(mem_fail_handler fail_handler) {
    fail_handler();
}

void free_mem_mg() {
    if (ctx.size == 0)
        return;
    m_reset_state(0);
    mem_free(ctx.ptrs);
    ctx.size = 0;
    ctx.idx = 0;
    ctx.ptrs = NULL;
}

enum mem_mg_init_res init_mem_mg(mem_fail_handler fail_handler) {
    if (ctx.size > 0)
        return ALREADY_INITIALIZED;
    intptr_t **ptrs = mem_alloc(DEFAULT_START_SIZE * sizeof(intptr_t));
    if (!ptrs) {
        handle_fail(fail_handler);
        return FAILED;
    }
    ctx.fail_handler = fail_handler;
    ctx.ptrs = ptrs;
    ctx.size = DEFAULT_START_SIZE;
    ctx.idx = 0;
    return OK;
}

static bool mem_mg_grow() {
    size_t nsize = ctx.size * 2;
    intptr_t **ptrs = mem_realloc((memptr_t) ctx.ptrs, sizeof(intptr_t) * nsize);
    if (!ptrs)
        return false;
    ctx.ptrs = ptrs;
    ctx.size = nsize;
    return true;
}

memptr_t m_alloc(size_t size) {
    if (ctx.idx == ctx.size && !mem_mg_grow())
        goto fail;
    memptr_t ptr = mem_alloc(size);
    if (!ptr)
        goto fail;
    ctx.ptrs[ctx.idx++] = ptr;
    return ptr;
fail:
    ctx.fail_handler();
    return NULL;
}

m_state m_current_state() {
    return ctx.idx;
}

void m_reset_state(m_state idx) {
    while (ctx.idx-- > idx)
        mem_free(ctx.ptrs[ctx.idx]);
    ctx.idx = idx;
}

void m_clear_states(m_state a, m_state b) {
    m_state tmp;
    if (a > b) {
        tmp = a;
        a = b;
        b = tmp;
    }

    for (m_state i = a; i < b; i++) {
        mem_free(ctx.ptrs[i]);
        ctx.ptrs[i] = NULL;
    }        
}
