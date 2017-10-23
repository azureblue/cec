#include <stdbool.h>
#include "cec_starter_split.h"
#include "cec_context.h"
#include "cec_starter_omp.h"

cec_mat ** group_points_by_clusters(cec_mat *x, int k, vec_i *clustering_vector) {
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

static cec_mat * create_sub_matrix(const cec_mat *src, int start_row, int end_row) {
    cec_mat * mat = cec_matrix_create(end_row - start_row, src->n);
    for (int i = start_row; i < end_row; i++)
        mat_copy_row(src, i, mat, i - start_row);
    return mat;
}

static cec_mat * create_single_row_matrix(int n, const double *data) {
    cec_mat *mat = cec_matrix_create(1, n);
    for (int i = 0; i < n; ++i)
        cec_matrix_set_element(mat, 0, i, data[i]);
    return mat;
}

double calculate_single_cluster_energy(const cec_mat *x_mat, const double* cluster_pos, int min_card, cec_model_spec model_spec) {
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

bool try_split_cluster(cec_mat *x_mat, double *cluster_pos, int min_card, int max_iters,
                              cec_model_spec model_spec, cec_mat *out_c_mat) {
    int n = x_mat->n;
    double single_cluster_energy = calculate_single_cluster_energy(x_mat,cluster_pos, min_card, model_spec);

    mem_state_id split_start = mem_track_start();
    vec_i *vc_two_centers = create_one_vc(2);

    cec_centers_par cen_double_par = {
            .centers_mat = NULL,
            .init_m = KMEANSPP,
            .var_centers = vc_two_centers
    };

    cec_mat *c_mat = cec_matrix_create(2, n);
    centers_init *split_ci = create_centers_init(&cen_double_par, x_mat->m);

    cec_control_par split_con_par = {
            .max_iterations = max_iters,
            .min_card = min_card,
            .starts = 10,
            .threads = 0
    };

    cec_model_spec two_specs[] = {model_spec, model_spec};
    cec_models_par split_models = {.len = 2, .model_specs = two_specs};

    cec_out *split_results;
    res_code split_res_code = cec_perform_starts(
            x_mat,
            c_mat,
            split_ci,
            &split_con_par,
            &split_models,
            &split_results
    );

    bool split_success = false;

    if (split_res_code == NO_ERROR
            && cec_final_centers_number(split_results) == 2
            && cec_final_energy(split_results) < single_cluster_energy) {
        cec_matrix_copy_data(split_results->centers, out_c_mat);
        split_success = true;
    }
    mem_reset_state(split_start);
    return split_success;
}

res_code run_single_start_with_given_clusters(cec_mat *x_mat, cec_mat *c_mat, int min_card, int max_iters,
                                                     cec_model_spec model_spec, cec_out ** results) {
    mem_state_id run_start = mem_track_start();
    vec_i *vc = create_one_vc(c_mat->m);

    cec_centers_par centers_final_par = {
            .centers_mat = c_mat,
            .init_m = NONE,
            .var_centers = vc
    };

    centers_init *ci = create_centers_init(&centers_final_par, x_mat->m);

    cec_control_par con_par = {
            .max_iterations = max_iters,
            .min_card = min_card,
            .starts = 1,
            .threads = 1
    };

    cec_model_spec model_specs_final[c_mat->m];
    cec_models_par mod_par = {.len = c_mat->m, .model_specs = model_specs_final};
    for (int i = 0; i < c_mat->m; ++i)
        model_specs_final[i] = model_spec;

    mem_state_range split_level_range = mem_track_end(run_start);
    res_code res = cec_perform_starts(x_mat,
                                      c_mat,
                                      ci,
                                      &con_par,
                                      &mod_par,
                                      results);
    mem_free_range(split_level_range);
    return res;
}

res_code cec_perform_split(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                        cec_out **results) {
    cec_model_spec model_spec = models->model_specs[0];
    int n = x_mat->n;
    int k = centers->var_centers->ar[0];

    int split_depth = 2;

    res_code current_res_code;
    mem_state_range current_res_range = mem_empty_range();
    cec_out *current_res;
    mem_track(
            current_res_code = cec_perform_starts(
                    x_mat,
                    cec_matrix_create(k, n),
                    create_centers_init(centers, x_mat->m),
                    control,
                    models,
                    &current_res),
            &current_res_range
    );

    if (current_res_code != NO_ERROR) {
        mem_free_range(current_res_range);
        return current_res_code;
    }

    for (int split_level = 0; split_level < split_depth; split_level++) {
        mem_state_id split_level_start = mem_track_start();
        k = current_res->centers->m;
        cec_mat **split = group_points_by_clusters(x_mat, k, current_res->clustering_vector);
        cec_mat *split_centers = cec_matrix_create(k * 2, n);
        cec_mat *split_res_c_mat = cec_matrix_create(2, n);
        int split_centers_k = k;

        for (int i = 0; i < k; i++) {
            if (split[i]->m == 0)
                continue;
            mat_copy_row(current_res->centers, i, split_centers, i);
            cec_mat *split_x_mat = split[i];
            bool split_success = try_split_cluster(split_x_mat,
                              mat_row(current_res->centers, i),
                              control->min_card,
                              control->max_iterations,
                              model_spec,
                              split_res_c_mat
            );
            if (split_success) {
                mat_copy_row(split_res_c_mat, 0, split_centers, i);
                mat_copy_row(split_res_c_mat, 1, split_centers, split_centers_k++);
            } else
                mat_copy_row(current_res->centers, i, split_centers, i);
        }

        cec_mat *level_c_mat = create_sub_matrix(split_centers, 0, split_centers_k);

        mem_state_range split_level_range = mem_track_end(split_level_start);
        mem_free_range(current_res_range);
        mem_track(run_single_start_with_given_clusters(x_mat, level_c_mat, control->min_card, control->max_iterations,
                          model_spec, &current_res),
                  &current_res_range);

        mem_free_range(split_level_range);
    }

    *results = current_res;

    return current_res_code;
}
