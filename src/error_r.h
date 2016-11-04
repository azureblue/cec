#ifndef ERROR_R_H
#define ERROR_R_H

#include "errors.h"
#if __STDC_VERSION__ >= 201112L
    #include <stdnoreturn.h>
#else
    #include <R_ext/Error.h>
    #define noreturn NORET
#endif

void noreturn error_r(enum error_code);

#endif /* ERROR_R_H */
