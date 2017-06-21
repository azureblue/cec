#include "init_utils_r.h"
#include "error_r.h"
#include "mem_mg.h"
#include "rand.h"

void noreturn error_r_mem_error() {
    cec_clean_env();
    error_r(MEM_ALLOC_ERROR);
};

void release_cec_mem_r() {
    free_mem_mg();
}

void cec_init_env() {
    init_mem_mg(error_r_mem_error);
    cec_rand_init();
}

void cec_clean_env() {
    cec_rand_end();
    free_mem_mg();
}
