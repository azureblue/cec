#include "cec_r.h"
#include "error_r.h"
#include "cec_params.h"
#include "cec_context.h"
#include "cec_params_r.h"
#include "result_r.h"
#include "cec_r_utils.h"
#include "init_utils_r.h"
#include "cec_starter_omp.h"

SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r)
{
    cec_init_env();
    int n = Rf_ncols(x);
    cec_centers_par *centers = get_centers_param(centers_param_r);
    cec_control_par *control = get_control_param(control_param_r);
    cec_models_par *models = get_models_param(models_param_r, n);
    cec_mat * x_mat = create_from_R_matrix(x);
    cec_out * out = create_cec_out_for_all_starts(x_mat, centers, control);
    res_code all_res = cec_perform(create_from_R_matrix(x), centers, control, models, out);

    SEXP result = NULL;

    if (all_res == NO_ERROR) {
        PROTECT(result = create_R_result(out));
        UNPROTECT(1);
        cec_clean_env();
        return result;
    }
    cec_clean_env();
    error_r(all_res);
    return NULL;
}
