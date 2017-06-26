#ifndef CEC_CEC_PARAMS_R_H
#define CEC_CEC_PARAMS_R_H

#include <Rdefines.h>
#include "cec_params.h"

cec_centers_par * get_centers_param(SEXP centers_param_r);
cec_control_par * get_control_param(SEXP control_param_r);
cec_models_par * get_models_param(SEXP models_param_r, int n);

#endif //CEC_CEC_PARAMS_R_H
