#include "alloc_utils_r.h"
#include "error_r.h"
#include "mem_mg.h"

void noreturn error_r_mem_error() {
    free_mem_mg();
    error_r(MEM_ALLOC_ERROR);
};

void release_cec_mem_r() {
    free_mem_mg();
}
