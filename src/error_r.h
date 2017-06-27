#ifndef ERROR_R_H
#define ERROR_R_H

#include "errors.h"
#include "noret.h"

void noreturn error_r(enum cec_result_code);
void noreturn missing_param_error_r(const char * param_name);
#endif /* ERROR_R_H */
