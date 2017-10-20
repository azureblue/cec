#include <omp.h>
#include "cec_starter_omp.h"
#include "cec.h"
#include "models/models.h"
#include "centers_init.h"

cec_out * create_cec_out_for_all_starts(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control) {
    int m = x_mat->m;
    int n = x_mat->n;
    int vc_max_k = max_i(centers->var_centers);
    return create_cec_output(m, vc_max_k, n, control->max_iterations);
}

static cec_mat * create_copy_from(cec_mat *c_mat, int k_rows) {
    cec_mat *res = cec_matrix_create(k_rows, c_mat->n);
    for (int i = 0; i < k_rows; ++i)
        cec_matrix_copy_row(c_mat, i, res, i);
    return res;
}

res_code cec_perform(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                     cec_out *results) {
    int vc_len = centers->var_centers->len;
    int m = x_mat->m;
    int n = x_mat->n;
    int vc_max_k = max_i(centers->var_centers);
    mem_state_id ms_start = mem_track_start();
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
    centers_init * ci = create_centers_init(centers, m);
    struct cec_model *cec_models[threads][models->len];
    cec_mat *local_c_mat[threads];
    cec_ctx *ctx[threads];

    for (int t = 0; t < threads; t++) {
        for (int i = 0; i < models->len; i++)
            cec_models[t][i] = create_model(&models->model_specs[i]);
    }

    for (int vc = 0; vc < vc_len; vc++) {
        int vc_k = centers->var_centers->ar[vc];
        c_mat->m = vc_k;
        mem_state_id vc_start = mem_track_start();
        for (int t = 0; t < threads; t++) {
            local_c_mat[t] = create_copy_from(c_mat, vc_k);
            ctx[t] = create_cec_context(
                    create_cec_input(x_mat, local_c_mat[t], cec_models[t], control->max_iterations, control->min_card),
                    create_cec_output(m, vc_max_k, n, control->max_iterations)
            );
        }

        #pragma omp parallel default(shared)
        {
            const int th = omp_get_thread_num();
            #pragma omp for nowait
            for (int start = 0; start < starts; start++) {

                #pragma omp critical
                cec_init_centers(ci, x_mat, local_c_mat[th]);

                res_code res = cec_start(ctx[th]);
                if (res)
                    continue;

                double energy = cec_final_energy(ctx[th]->results);
                // best_energy may be stale but never lower than the true value
                if (energy >= best_energy)
                    continue;
                #pragma omp critical
                    // need to check again!
                if (energy < best_energy) {
                    cec_copy_results_content(ctx[th]->results, results);
                    best_energy = energy;
                    all_res = NO_ERROR;
                }
            }
        }
        mem_reset_state(vc_start);
    }
    omp_set_num_threads(threads_default);
    mem_reset_state(ms_start);
    return all_res;
}
