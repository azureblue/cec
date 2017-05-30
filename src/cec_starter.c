#include <omp.h>
#include "errors.h"
#include "cec_context.h"
#include "cec_r.h"
#include "cec.h"

struct cec_model * create_model(struct cec_model_spec * model_spec)
{
    struct cec_model *model = alloc(struct cec_model);
    switch (model_spec->type)
    {
        case ALL:
            model->cross_entropy = h_all;
            model->cross_entropy_context = create_cross_entropy_context_all(model_spec->n);
            break;
        case SPHERICAL:
            model->cross_entropy = h_spherical;
            model->cross_entropy_context = create_cross_entropy_context_spherical();
            break;
        case DIAGONAL:
            model->cross_entropy = h_diagonal;
            model->cross_entropy_context = create_cross_entropy_context_diagonal();
            break;
        case FIXED_R:
            model->cross_entropy = h_fixed_r;
            model->cross_entropy_context = create_cross_entropy_context_fixedr(model_r_params(model_spec)->r);
            break;
        case GIVEN_COVARIANCE:
            model->cross_entropy = h_given_covariance;
            model->cross_entropy_context = create_cross_entropy_context_covariance(
                    model_covariances_params(model_spec)->cov, model_covariances_params(model_spec)->cov_inv);
            break;
        case FIXEDEIGENVALUES:
            model->cross_entropy = h_fixedeigenvalues;
            model->cross_entropy_context = create_cross_entropy_context_eigenvalues(
                    model_eigenvalues_params(model_spec)->given_eigenvalues->len,
                    model_eigenvalues_params(model_spec)->given_eigenvalues->ar);
    }

    return model;
}


res_code cec_perform(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                     cec_out **results) {
    int vc_len = centers->var_centers->len;
    int m = x_mat->m;
    int n = x_mat->n;
    int vc_max_k = -1;
    for (int i = 0; i < vc_len; i++)
        vc_max_k = vc_max_k < centers->var_centers->ar[i] ? centers->var_centers->ar[i] : vc_max_k;
    cec_out *best_result = create_cec_output(m, vc_max_k, n, control->max_iterations);
    m_state ms_start = m_current_state();
    cec_mat * c_mat = cec_matrix_create(vc_max_k, n);
    if (centers->centers_mat)
        cec_matrix_copy_data(centers->centers_mat, c_mat);
    double best_energy = BIG_DOUBLE;
    res_code all_res = UNKNOWN_ERROR;
    int starts = control->starts;

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
    m_reset_state(ms_start);
    *results = best_result;
    return all_res;
}

