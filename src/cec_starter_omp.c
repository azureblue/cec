#include <omp.h>
#include "cec_starter_omp.h"
#include "cec.h"
#include "models/models.h"

static int max_i(const vec_i *ar) {
    int len = ar->len;
    int max = INT_MIN;
    for (int i = 0; i < len; i++)
        max = max > ar->ar[i] ? max : ar->ar[i];
    return max;
}

res_code cec_perform(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                     cec_out **results) {
    int vc_len = centers->var_centers->len;
    int m = x_mat->m;
    int n = x_mat->n;
    int vc_max_k = max_i(centers->var_centers);
    cec_out *best_result = create_cec_output(m, vc_max_k, n, control->max_iterations);
    m_state ms_start = m_current_state();
    cec_mat * c_mat = cec_matrix_create(vc_max_k, n);
    if (centers->centers_mat)
        cec_matrix_copy_data(centers->centers_mat, c_mat);
    double best_energy = BIG_DOUBLE;
    res_code all_res = UNKNOWN_ERROR;
    int starts = control->starts;
    int threads_default = omp_get_max_threads();
    int threads = control->threads == 0 ? threads_default : control->threads;

    if (threads > starts)
        threads = starts;

    omp_set_num_threads(threads);

    for (int vc = 0; vc < vc_len; vc++) {
        int vc_k = centers->var_centers->ar[vc];
        c_mat->m = vc_k;
        m_state vc_start = m_current_state();
#pragma omp parallel default(shared)
        {
            cec_mat *local_c_mat;
            cec_in *in;
            cec_out *out;
            cec_ctx *ctx;
            struct cec_model *cec_models[models->len];
#pragma omp critical
            {
                for (int i = 0; i < models->len; i++)
                    cec_models[i] = create_model(&models->model_specs[i]);

                local_c_mat = cec_matrix_create_copy(c_mat);
                in = create_cec_input(x_mat, local_c_mat, cec_models, control->max_iterations, control->min_card);
                out = create_cec_output(m, vc_k, n, control->max_iterations);
                ctx = create_cec_context(in, out);
            }
#pragma omp for nowait
            for (int start = 0; start < starts; start++) {

#pragma omp critical
                cec_init_centers(x_mat, local_c_mat, centers->init_m);

                int res = cec_start(ctx);

                if (res == NO_ERROR) {
                    double energy = cec_final_energy(out);
                    // best_energy may be stale but never lower than the true value
                    if (energy < best_energy) {
#pragma omp critical
                        // need to check again!
                        if (energy < best_energy) {
                            cec_copy_results_content(out, best_result);
                            best_energy = energy;
                            all_res = NO_ERROR;
                        }
                    }
                }
            }
        }
        m_reset_state(vc_start);
    }
    omp_set_num_threads(threads_default);
    m_reset_state(ms_start);
    *results = best_result;
    return all_res;
}

