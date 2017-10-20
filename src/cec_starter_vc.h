#ifndef CEC_CEC_STARTER_VC_H
#define CEC_CEC_STARTER_VC_H

#include "errors.h"
#include "cec_params.h"
#include "cec_context.h"

res_code cec_perform_vc(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                        cec_out **results);

#endif //CEC_CEC_STARTER_VC_H
