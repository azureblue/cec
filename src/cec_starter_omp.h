#ifndef CEC_CEC_STARTER_OMP_H
#define CEC_CEC_STARTER_OMP_H

#include "errors.h"
#include "matrix.h"
#include "cec_params.h"
#include "cec_context.h"

res_code cec_perform(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control,
                     cec_models_par *models, cec_out **results);

#endif //CEC_CEC_STARTER_OMP_H
