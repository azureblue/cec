#include <omp.h>
#include "cec_starter_omp.h"
#include "cec.h"
#include "models/models.h"

typedef struct {
    cec_in * in;
    cec_out * out;
    cec_tmp_data * temp_data;
} cec_ctx;

res_code cec_perform_starts(const cec_mat *x_mat, cec_mat *c_mat, centers_initializer *ci, initial_assignment *ia,
                            cec_control_par *control,
                            cec_models_par *models,
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

    const cec_model *local_models[threads][models->len];
    cec_mat *local_centers[threads];
    vec_i *local_assignments[threads];
    cec_ctx *ctx[threads];

    for (int t = 0; t < threads; t++) {
        ctx[t] = alloc(cec_ctx);
        *(ctx[t]->in = alloc(cec_in)) = (cec_in) {
                .points = x_mat,
                .centers = local_centers[t] = cec_matrix_create(k, n),
                .initial_assignment = local_assignments[t] = vec_i_create(m),
                .models = local_models[t],
                .max_iterations = control->max_iterations,
                .min_card = control->min_card
        };
        ctx[t]->out = create_cec_output(m, k, n, control->max_iterations);
        ctx[t]->temp_data = create_temp_data(k, n);

        for (int i = 0; i < models->len; i++)
            local_models[t][i] = create_model(&models->model_specs[i]);
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
            cec_init_centers(ci, x_mat, local_centers[th]);

            cec_assign_points(ia, x_mat, local_centers[th], local_assignments[th]);

            if(cec_start(ctx[th]->in, ctx[th]->out, ctx[th]->temp_data) != NO_ERROR)
                continue;

            double energy = cec_final_energy(ctx[th]->out);
            // best_energy may be stale but never lower than the true value
            if (energy >= best_energy)
                continue;
            #pragma omp critical
            // need to check again!
            if (energy < best_energy) {
                cec_copy_results_content(ctx[th]->out, best);
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
