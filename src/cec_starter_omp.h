#ifndef CEC_STARTER_OMP_H
#define CEC_STARTER_OMP_H

#include "errors.h"
#include "matrix.h"
#include "cec_params.h"
#include "cec_context.h"
#include "centers_init.h"

res_code cec_perform(cec_mat *x_mat, cec_mat *c_mat, centers_init *ci, cec_control_par *control, cec_models_par *models,
                     cec_out **results);

#endif //CEC_STARTER_OMP_H
