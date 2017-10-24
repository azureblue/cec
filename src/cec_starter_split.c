#include <stdbool.h>
#include <stdlib.h>
#include "cec_starter_split.h"
#include "cec_starter_omp.h"
#include "cec_r_utils.h"
#include "kmeanspp.h"

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

void debug_clustering(cec_mat *x_mat, int *clustering) {
    FILE *file = fopen("/media/docs_ram/debug.R", "w");
    fprintf(file, "X11(); \n plot(matrix(c(");

    for (int j = 0; j < x_mat->n - 1; ++j)
        for (int i = 0; i < x_mat->m; ++i)
            fprintf(file, "%lf, ", cec_matrix_element(x_mat, i, j));

    for (int j = 0; j < x_mat->m - 1; ++j)
        fprintf(file, "%lf, ", cec_matrix_element(x_mat, j, x_mat->n - 1));

    fprintf(file, "%lf", cec_matrix_element(x_mat, x_mat->m - 1, x_mat->n - 1));

    fprintf(file, "),,2), col=c(");
    for (int i = 0; i < x_mat->m - 1; ++i) {
        fprintf(file, "%d, ", clustering[i] + 1);
    }
    fprintf(file, "%d", clustering[x_mat->m - 1] + 1);
    fprintf(file, "), xlab=1, ylab=2, pch=20); \nSys.sleep(2);\n while (!is.null(dev.list())) Sys.sleep(0.2);");
    fclose(file);

    system("Rscript /media/docs_ram/debug.R");
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

bool try_split_cluster(cec_mat *x_mat, double *cluster_pos, int min_card, int max_iters,
                       cec_model_spec model_spec, cec_mat *out_c_mat, vec_i *assignment) {
    int n = x_mat->n;
    double single_cluster_energy = calculate_single_cluster_energy(x_mat, cluster_pos, min_card, model_spec);

    mem_state_id split_start = mem_track_start();
    cec_mat *c_mat = cec_matrix_create(2, n);

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

    mem_state_range split_level_range = mem_track_end(run_start);

    res_code res = cec_perform_starts(x_mat,
                                      c_mat,
                                      create_copy_initializer(c_mat),
                                      create_initial_assignment_fixed(initial_cluster),
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

    int split_depth = 4;

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

    for (int split_level = 0; split_level < split_depth; split_level++) {
        vec_i_copy(current_res->clustering_vector, cluster);
        mem_state_id split_level_start = mem_track_start();
        k = current_res->centers->m;
        cluster_split *split = group_points_by_clusters(x_mat, k, current_res->clustering_vector);
        cec_mat *split_centers = cec_matrix_create(k * 2, n);
        cec_mat *split_res_c_mat = cec_matrix_create(2, n);
        int split_centers_k = k;

        for (int i = 0; i < k; i++) {
            int split_m = split[i].points->m;
            if (split_m == 0)
                continue;
            mat_copy_row(current_res->centers, i, split_centers, i);
            cec_mat *split_x_mat = split[i].points;
            vec_i *mapping = split[i].mapping;
            vec_i *split_cluster = vec_i_create(split_x_mat->m);
            bool split_success = try_split_cluster(split_x_mat,
                                                   mat_row(current_res->centers, i),
                                                   control->min_card,
                                                   control->max_iterations,
                                                   model_spec,
                                                   split_res_c_mat,
                                                   split_cluster
            );
            if (split_success) {
                mat_copy_row(split_res_c_mat, 0, split_centers, i);
                mat_copy_row(split_res_c_mat, 1, split_centers, split_centers_k);
                for (int p = 0; p < split_m; p++) {
                    int c = split_cluster->ar[p];
                    if (c == 0)
                        cluster->ar[mapping->ar[p]] = i;
                    else
                        cluster->ar[mapping->ar[p]] = split_centers_k;
                }
                split_centers_k++;

            } else
                mat_copy_row(current_res->centers, i, split_centers, i);
        }
        //debug_clustering(x_mat, cluster->ar);
        cec_mat *level_c_mat = create_sub_matrix(split_centers, 0, split_centers_k);

        mem_state_range split_level_range = mem_track_end(split_level_start);

        mem_free_range(current_res_range);
        res_code rc;
        mem_track(&current_res_range,
                  rc = run_single_start_with_given_clusters(x_mat, level_c_mat, cluster, control->min_card,
                                                       control->max_iterations, model_spec, &current_res));

        printf("res: %d, iters: %d\n", rc, current_res->iterations);
        debug_clustering(x_mat, current_res->clustering_vector->ar);

        mem_free_range(split_level_range);
    }

    *results = current_res;

    return current_res_code;
}
