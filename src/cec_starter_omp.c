#include <omp.h>
#include "cec_starter_omp.h"
#include "cec.h"
#include "models/models.h"
#include "init_utils_r.h"
#include "centers_init.h"

res_code cec_perform(cec_mat *x_mat, cec_mat *c_mat, centers_init *ci, cec_control_par *control, cec_models_par *models,
                     cec_out **results) {
    mem_state_id ms_start = mem_track_start();
    int m = x_mat->m;
    int k = c_mat->m;
    int n = x_mat->n;
    int starts = control->starts;
    int threads_default = omp_get_max_threads();
    int threads = control->threads;
    if (threads == 0)
        threads = threads_default;
    if (threads > starts)
        threads = starts;

    omp_set_num_threads(threads);

    struct cec_model *cec_models[threads][models->len];
    cec_mat *local_c_mat[threads];
    cec_ctx *ctx[threads];

    for (int t = 0; t < threads; t++) {
        local_c_mat[t] = cec_matrix_create(k, n);
        ctx[t] = create_cec_context(
                create_cec_input(x_mat, local_c_mat[t], cec_models[t], control->max_iterations, control->min_card),
                create_cec_output(m, k, n, control->max_iterations)
        );

        for (int i = 0; i < models->len; i++)
            cec_models[t][i] = create_model(&models->model_specs[i]);
    }

    cec_out * best = create_cec_output(m, k, n, control->max_iterations);
    double best_energy = BIG_DOUBLE;
    res_code all_res = UNKNOWN_ERROR;

    #pragma omp parallel default(shared)
    {
        const int th = omp_get_thread_num();
        #pragma omp for nowait
        for (int start = 0; start < starts; start++) {
            #pragma omp critical
            cec_init_centers(ci, x_mat, local_c_mat[th]);

            if(cec_start(ctx[th]) != NO_ERROR)
                continue;

            double energy = cec_final_energy(ctx[th]->results);
            // best_energy may be stale but never lower than the true value
            if (energy >= best_energy)
                continue;
            #pragma omp critical
            // need to check again!
            if (energy < best_energy) {
                cec_copy_results_content(ctx[th]->results, best);
                best_energy = energy;
                all_res = NO_ERROR;
            }
        }
    }

    omp_set_num_threads(threads_default);
    mem_state_range ms_range = mem_track_end(ms_start);

    if (all_res == NO_ERROR) {
        *results = create_cec_output(m, best->centers->m, n, control->max_iterations);
        cec_copy_results_content(best, *results);
    }
    mem_free_range(ms_range);
    return all_res;
}
