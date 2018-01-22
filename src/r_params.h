#ifndef CEC_PARAMS_R_H
#define CEC_PARAMS_R_H

#include "params.h"
#include <Rdefines.h>

namespace cec {
    centers_param get_centers_param(SEXP centers_param_r);

    control_param get_control_param(SEXP control_param_r);

    models_param get_models_param(SEXP models_param_r, int n);

    split_param get_split_param(SEXP split_param_r);
}

#endif //CEC_PARAMS_R_H
