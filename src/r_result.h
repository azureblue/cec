#ifndef RESULT_R_H
#define RESULT_R_H

#include <Rdefines.h>
#include "starter.h"

namespace cec {
    SEXP create_R_result(const clustering_results &res);
}

#endif //RESULT_R_H
