#ifndef RESULT_R_H
#define RESULT_R_H

#include "starter.h"

#include <Rdefines.h>

namespace cec {
    SEXP create_R_result(const clustering_results &res);
}

#endif //RESULT_R_H
