#ifndef CEC_PARAMS_R_H
#define CEC_PARAMS_R_H

#include "cec_params.h"

extern "C" {
#include <Rdefines.h>
};

namespace cec {
    cec_centers_param get_centers_param(SEXP centers_param_r);

    cec_control_param get_control_param(SEXP control_param_r);

    cec_models_param get_models_param(SEXP models_param_r, int n);
}

#endif //CEC_PARAMS_R_H
