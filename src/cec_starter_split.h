#include "errors.h"
#include "cec_params.h"
#include "cec_context.h"

#ifndef CEC_SPLIT_H
#define CEC_SPLIT_H

res_code cec_perform_split(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                           cec_out **results);

#endif //CEC_SPLIT_H
