#include <stdbool.h>
#include <stdlib.h>
#include "cec_starter_split.h"
#include "cec_starter_omp.h"
#include "cec_r_utils.h"
#include "kmeanspp.h"
#include "cec.h"

typedef struct {
    cec_mat *points;
    vec_i *mapping;
} cluster_split;

cluster_split *group_points_by_clusters(cec_mat *x, int k, vec_i *clustering_vector) {
    int m = clustering_vector->len;
    int n = x->n;
    int cluster_sizes[k];
    int cluster_idxs[k];
    cluster_split *split = alloc_n(cluster_split, k);
    for (int i = 0; i < k; i++) {
        cluster_sizes[i] = 0;
        cluster_idxs[i] = 0;
    }
    for (int i = 0; i < m; i++)
        cluster_sizes[clustering_vector->ar[i]]++;
    for (int i = 0; i < k; i++) {
        split[i].points = mat_create(cluster_sizes[i], n);
        split[i].mapping = vec_i_create(cluster_sizes[i]);
    }
    for (int i = 0; i < m; i++) {
        int dest_cluster = clustering_vector->ar[i];
        int idx = cluster_idxs[dest_cluster];
        split[dest_cluster].mapping->ar[idx] = i;
        mat_copy_row(x, i, split[dest_cluster].points, idx);
        cluster_idxs[dest_cluster]++;
    }
    return split;
}

bool same_pos(const double *a, const double *b, int n) {
    for (int i = 0; i < n; i++)
        if (a[i] != b[i])
            return false;
    return true;
}

static vec_i *create_one_vc(int k) {
    vec_i *vc = vec_i_create(1);
    vc->ar[0] = k;
    return vc;
}

static cec_mat *create_sub_matrix(const cec_mat *src, int start_row, int end_row) {
    cec_mat *mat = cec_matrix_create(end_row - start_row, src->n);
    for (int i = start_row; i < end_row; i++)
        mat_copy_row(src, i, mat, i - start_row);
    return mat;
}

static cec_mat *create_single_row_matrix(int n, const double *data) {
    cec_mat *mat = cec_matrix_create(1, n);
    for (int i = 0; i < n; ++i)
        cec_matrix_set_element(mat, 0, i, data[i]);
    return mat;
}

double calculate_single_cluster_energy(const cec_mat *x_mat, const double *cluster_pos, int min_card,
                                       cec_model_spec model_spec) {
    mem_state_id start = mem_track_start();
    cec_mat *c_mat = create_single_row_matrix(x_mat->n, cluster_pos);
    vec_i *single_vc = create_one_vc(1);
    cec_centers_par cen_par = {
            .centers_mat = c_mat,
            .init_m = NONE,
            .var_centers = single_vc
    };

    centers_initializer *ci = create_centers_init(&cen_par, x_mat->m);
    initial_assignment *ia = create_initial_assignment_closest();

    cec_control_par con_par = {
            .max_iterations = 1,
            .min_card = min_card,
            .starts = 1,
            .threads = 1
    };

    cec_model_spec model_specs[] = {model_spec};
    cec_models_par mod_par = {.len = 1, .model_specs = model_specs};

    cec_out *results;
    cec_perform_starts(x_mat, c_mat, ci, ia, &con_par, &mod_par, &results);

    double energy = cec_final_energy(results);
    mem_reset_state(start);
    return energy;
}

bool try_split_cluster(cec_mat *x_mat, double *cluster_pos, int min_card, int max_iters, int tries,
                       cec_model_spec model_spec, cec_mat *out_c_mat, vec_i *assignment) {
    int n = x_mat->n;
    double single_cluster_energy = calculate_single_cluster_energy(x_mat, cluster_pos, min_card, model_spec);

    mem_state_id split_start = mem_track_start();
    cec_mat *c_mat = cec_matrix_create(2, n);

    cec_control_par split_con_par = {
            .max_iterations = max_iters,
            .min_card = min_card,
            .starts = tries,
            .threads = 0
    };

    cec_model_spec two_specs[] = {model_spec, model_spec};
    cec_models_par split_models = {.len = 2, .model_specs = two_specs};

    cec_out *split_results;
    res_code split_res_code = cec_perform_starts(
            x_mat,
            c_mat,
            create_kmeanspp_initializer(x_mat->m),
            create_initial_assignment_closest(),
            &split_con_par,
            &split_models,
            &split_results
    );

    bool split_success = false;
    //debug_clustering(x_mat, split_results->clustering_vector->ar);
    if (split_res_code == NO_ERROR
        && cec_final_centers_number(split_results) == 2
        && cec_final_energy(split_results) < single_cluster_energy) {
        mat_copy_data(split_results->centers, out_c_mat);
        vec_i_copy(split_results->clustering_vector, assignment);

        split_success = true;
    }
    mem_reset_state(split_start);
    return split_success;
}

res_code run_single_start_with_given_clusters(cec_mat *x_mat, cec_mat *c_mat, vec_i *initial_cluster, int min_card,
                                              int max_iters, cec_model_spec model_spec, cec_out **results) {
    mem_state_id run_start = mem_track_start();

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

    mem_state_range run_range = mem_track_end(run_start);

    res_code res = cec_perform_starts(x_mat,
                                      c_mat,
                                      create_copy_initializer(c_mat),
                                      create_initial_assignment_fixed(initial_cluster),
                                      &con_par,
                                      &mod_par,
                                      results);
    mem_free_range(run_range);
    return res;
}

static res_code cec_perform_split_start(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                                        cec_split_par *split_par, cec_out **results) {
    cec_model_spec model_spec = models->model_specs[0];
    int n = x_mat->n;
    int k = centers->var_centers->ar[0];

    int split_depth = split_par->depth;

    res_code current_res_code;
    mem_state_range current_res_range = mem_empty_range();
    vec_i *cluster = vec_i_create(x_mat->m);
    cec_out *current_res;
    mem_track(
            &current_res_range,
            current_res_code = cec_perform_starts(
                    x_mat,
                    cec_matrix_create(k, n),
                    create_centers_init(centers, x_mat->m),
                    create_initial_assignment_closest(),
                    control,
                    models,
                    &current_res)
    );

    if (current_res_code != NO_ERROR) {
        mem_free_range(current_res_range);
        return current_res_code;
    }

    mem_state_id moved_start = mem_track_start();
    bool *moved = alloc_n(bool, k);
    mem_state_range moved_range = mem_track_end(moved_start);
    for (int i = 0; i < k; i++)
        moved[i] = true;

    for (int split_level = 0; split_level < split_depth; split_level++) {
        vec_i_copy(current_res->clustering_vector, cluster);
        mem_state_id split_level_start = mem_track_start();
        k = current_res->centers->m;
        cluster_split *split = group_points_by_clusters(x_mat, k, current_res->clustering_vector);
        cec_mat *split_centers = mat_create(k * 2, n);
        cec_mat *split_res_c_mat = mat_create(2, n);
        int k_s = 0;
        bool splitted[k * 2];
        bool split_flag = false;

        for (int i = 0; i < k; i++) {
            int split_m = split[i].points->m;
            if (split_m == 0)
                continue;
            mat_copy_row(current_res->centers, i, split_centers, k_s);
            bool split_success = false;

            cec_mat *split_x_mat = split[i].points;
            vec_i *mapping = split[i].mapping;
            vec_i *split_assignment = vec_i_create(split_x_mat->m);
            if (moved[i])
                split_success = try_split_cluster(split_x_mat,
                                                   mat_row(current_res->centers, i),
                                                   control->min_card,
                                                   control->max_iterations,
                                                   split_par->tries,
                                                   model_spec,
                                                   split_res_c_mat,
                                                   split_assignment
            );
            split_flag = split_flag || split_success;
            if (split_success) {
                splitted[k_s] = true;
                splitted[k_s + 1] = true;
                mat_copy_row(split_res_c_mat, 0, split_centers, k_s);
                mat_copy_row(split_res_c_mat, 1, split_centers, k_s + 1);
                for (int p = 0; p < split_m; p++) {
                    if (split_assignment->ar[p] == 0)
                        cluster->ar[mapping->ar[p]] = k_s;
                    else
                        cluster->ar[mapping->ar[p]] = k_s + 1;
                }
                k_s += 2;

            } else {
                splitted[k_s] = false;
                mat_copy_row(current_res->centers, i, split_centers, k_s);
                for (int p = 0; p < split_m; p++)
                    cluster->ar[mapping->ar[p]] = k_s;
                k_s++;
            }
        }
        cec_mat *level_c_mat = create_sub_matrix(split_centers, 0, k_s);

        mem_state_range split_level_range = mem_track_end(split_level_start);

        res_code rc;
        bool need_another_split = false;
        if (split_flag) {
            mem_track(&current_res_range,
                      rc = run_single_start_with_given_clusters(x_mat, level_c_mat, cluster, control->min_card,
                                                                control->max_iterations, model_spec, &current_res));

            mem_track(&moved_range, moved = alloc_n(bool, k_s));

            for (int i = 0; i < x_mat->m; i++) {
                if (cluster->ar[i] != current_res->clustering_vector->ar[i]) {
                    moved[cluster->ar[i]] = true;
                    moved[current_res->clustering_vector->ar[i]] = true;
                }
            }
            for (int i = 0; i < k_s; i++) {
                if (splitted[i])
                    moved[i] = true;
                need_another_split = (need_another_split || moved[i]);
            }
        }

        mem_free_range(split_level_range);
        if (!need_another_split)
            break;
    }

    *results = current_res;

    return current_res_code;
}

res_code cec_perform_split(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *models,
                           cec_split_par *split, cec_out **results) {
    cec_out *best = NULL;
    res_code best_res_code = UNKNOWN_ERROR;
    double best_energy = BIG_DOUBLE;
    mem_state_range cec_out_best_mem_range = mem_empty_range();
    for (int i = 0; i < control->starts; i++){
        cec_out *split_start_res;
        mem_state_range local_res_range = mem_empty_range();
        res_code res;
        mem_track(&local_res_range,
                  res = cec_perform_split_start(x_mat, centers, control, models, split, &split_start_res));
        if (res == NO_ERROR && cec_final_energy(split_start_res) < best_energy) {
            mem_free_range(cec_out_best_mem_range);
            cec_out_best_mem_range = local_res_range;
            best = split_start_res;
            best_energy = cec_final_energy(split_start_res);
            best_res_code = res;
        } else
            mem_free_range(local_res_range);
    }
    *results = best;
    return best_res_code;
}
