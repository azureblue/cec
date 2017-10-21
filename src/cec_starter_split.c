#include <stdbool.h>
#include "cec_starter_split.h"
#include "cec_context.h"
#include "cec_starter_omp.h"

static cec_mat ** group_points_by_clusters(cec_mat *x, int k, vec_i *clustering_vector) {
    int m = clustering_vector->len;
    int n = x->n;
    int cluster_sizes[k];
    cec_mat **split = alloc_n(cec_mat*, k);
    for (int i = 0; i < k; i++)
        cluster_sizes[i] = 0;
    for (int i = 0; i < m; i++)
        cluster_sizes[clustering_vector->ar[i]]++;
    for (int i = 0; i < k; i++)
        split[i] = cec_matrix_create(cluster_sizes[i], n);
    for (int i = 0; i < m; i++) {
        int dest_cluster = clustering_vector->ar[i];
        cec_mat *dest_split = split[dest_cluster];
        array_copy(cec_matrix_const_row(x, i), cec_matrix_row(dest_split, --cluster_sizes[dest_cluster]), n);
    }
    return split;
}

static vec_i * create_one_vc(int k) {
    vec_i * vc = vec_i_create(1);
    vc->ar[0] = k;
    return vc;
}

static cec_mat * create_single_row_matrix(int n, double data[]) {
    cec_mat *mat = cec_matrix_create(1, n);
    for (int i = 0; i < n; ++i)
        cec_matrix_set_element(mat, 0, i, data[i]);
    return mat;
}

static double calculate_single_cluster_energy(cec_mat *x_mat, double* cluster_pos, int min_card, cec_model_spec model_spec) {
    mem_state_id start = mem_track_start();
    cec_mat *c_mat = create_single_row_matrix(x_mat->n, cluster_pos);
    vec_i * single_vc = create_one_vc(1);
    cec_centers_par cen_par = {
            .centers_mat = c_mat,
            .init_m = NONE,
            .var_centers = single_vc
    };

    centers_init *ci = create_centers_init(&cen_par, x_mat->m);

    cec_control_par con_par = {
            .max_iterations = 1,
            .min_card = min_card,
            .starts = 1,
            .threads = 1
    };

    cec_model_spec model_specs[] = {model_spec};
    cec_models_par mod_par = {.len = 1, .model_specs = model_specs};

    cec_out *results;
    cec_perform_starts(x_mat, c_mat, ci, &con_par, &mod_par, &results);

    double energy = cec_final_energy(results);
    mem_reset_state(start);
    return energy;
}

res_code cec_perform_split(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                        cec_out **results) {
    mem_state_id perform_split_start = mem_track_start();
    cec_model_spec model_spec = models->model_specs[0];
    int n = x_mat->n;
    int k = centers->var_centers->ar[0];

    centers_init *ci = create_centers_init(centers, x_mat->m);
    cec_out *initial_results;

    res_code initial_res = cec_perform_starts(x_mat, cec_matrix_create(k, n), ci, control, models, &initial_results);
    if (initial_res != NO_ERROR) {
        mem_reset_state(perform_split_start);
        return initial_res;
    }
    cec_mat ** split = group_points_by_clusters(x_mat, k, initial_results->clustering_vector);
    cec_mat * split_centers = cec_matrix_create(k * 2, n);
    int split_centers_k = k;

    for (int i = 0; i < k; i++) {
        if (split[i]->m == 0)
            continue;
        cec_matrix_copy_row(initial_results->centers, i, split_centers, i);
        cec_mat *split_x_mat = split[i];
        double single_cluster_energy
                = calculate_single_cluster_energy(split_x_mat, cec_matrix_row(initial_results->centers, i), control->min_card, model_spec);


        mem_state_id split_start = mem_track_start();
        vec_i * vc_two_centers = create_one_vc(2);

        cec_centers_par cen_double_par = {
                .centers_mat = NULL,
                .init_m = KMEANSPP,
                .var_centers = vc_two_centers
        };

        cec_mat *split_c_mat = cec_matrix_create(2, n);
        centers_init *split_ci = create_centers_init(&cen_double_par, split_x_mat->m);

        cec_control_par split_con_par = {
                .max_iterations = control->max_iterations,
                .min_card = control->min_card,
                .starts = 10,
                .threads = 0
        };

        cec_model_spec two_specs[] = {model_spec, model_spec};
        cec_models_par split_models = {.len = 2, .model_specs = two_specs};

        cec_out *split_results;
        res_code split_res_code = cec_perform_starts(
                split_x_mat,
                split_c_mat,
                split_ci,
                &split_con_par,
                &split_models,
                &split_results
        );

        if (split_res_code != NO_ERROR) {
            continue;
        }

        if (cec_final_centers_number(split_results) == 2 && cec_final_energy(split_results) < single_cluster_energy) {
            cec_matrix_copy_row(split_results->centers, 0, split_centers, i);
            cec_matrix_copy_row(split_results->centers, 1, split_centers, split_centers_k++);
        }
        mem_reset_state(split_start);
    }

    cec_model_spec model_specs_final[split_centers_k];
    cec_models_par mod_par_final = {.len = split_centers_k, .model_specs = model_specs_final};
    cec_mat *centers_final = cec_matrix_create(split_centers_k, n);
    for (int i = 0; i < split_centers_k; ++i) {
        model_specs_final[i] = model_spec;
        cec_matrix_copy_row(split_centers, i, centers_final, i);
    }
    vec_i * vc_final = create_one_vc(split_centers_k);

    cec_centers_par centers_final_par = {
            .centers_mat = centers_final,
            .init_m = NONE,
            .var_centers = vc_final
    };

    centers_init *ci_final = create_centers_init(&centers_final_par, x_mat->m);

    cec_control_par con_par_final = {
            .max_iterations = control->max_iterations,
            .min_card = control->min_card,
            .starts = 1,
            .threads = 1
    };

    mem_state_range split_range = mem_track_end(perform_split_start);

    res_code final_res = cec_perform_starts(
            x_mat,
            centers_final,
            ci_final,
            &con_par_final,
            &mod_par_final,
            results
    );

    mem_free_range(split_range);
    return final_res;
}
