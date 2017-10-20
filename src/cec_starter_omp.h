#ifndef CEC_STARTER_OMP_H
#define CEC_STARTER_OMP_H

#include "errors.h"
#include "matrix.h"
#include "cec_params.h"
#include "cec_context.h"

cec_out * create_cec_out_for_all_starts(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control);

res_code cec_perform(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control,
                     cec_models_par *models, cec_out *results);

#endif //CEC_STARTER_OMP_H
