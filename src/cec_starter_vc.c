#include "cec_starter_vc.h"
#include "cec_starter_omp.h"
#include "cec.h"

res_code cec_perform_vc(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                        cec_out **results) {
    mem_state_id vc_start = mem_track_start();
    int m = x_mat->m;
    int n = x_mat->n;
    cec_out *best = NULL;
    double best_energy = BIG_DOUBLE;
    mem_state_range cec_out_best_mem_range = mem_empty_range();
    centers_init *ci = create_centers_init(centers, x_mat->m);
    res_code all_res = UNKNOWN_ERROR;

    for (int vc = 0; vc < centers->var_centers->len; vc++) {
        int vc_k = centers->var_centers->ar[vc];
        mem_state_id vc_loop_start = mem_track_start();
        cec_mat *c_mat = cec_matrix_create(vc_k, n);
        cec_out *vc_out;

        if (cec_perform(x_mat, c_mat, ci, control, models, &vc_out) != NO_ERROR)
            continue;

        double energy = cec_final_energy(vc_out);
        mem_state_range vc_loop_range = mem_track_end(vc_loop_start);

        if (energy < best_energy) {
            all_res = NO_ERROR;
            mem_free_range(cec_out_best_mem_range);
            mem_track(best = create_cec_output(m, vc_k, n, control->max_iterations), &cec_out_best_mem_range);
            cec_copy_results_content(vc_out, best);
            best_energy = energy;
        }
        mem_free_range(vc_loop_range);
    }
    mem_state_range vc_range = mem_track_end(vc_start);

    if (all_res == NO_ERROR) {
        *results = create_cec_output(m, best->centers->m, n, control->max_iterations);
        cec_copy_results_content(best, *results);
    }
    mem_free_range(vc_range);
    return all_res;
}
