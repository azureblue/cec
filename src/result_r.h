#ifndef RESULT_R_H
#define RESULT_R_H

#include "single_start_input.h"

#include <Rdefines.h>


namespace cec {
    SEXP create_R_result(const single_start_results &res);
}

#endif //RESULT_R_H
