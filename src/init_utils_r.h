#ifndef ALLOC_UTILS_R_H
#define ALLOC_UTILS_R_H

#include "errors.h"
#include "noret.h"

void noreturn error_r_mem_error();
void release_cec_mem_r();
void cec_init_env();
void cec_clean_env();

#endif //ALLOC_UTILS_R_H
