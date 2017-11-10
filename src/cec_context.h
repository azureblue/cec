#ifndef CEC_CONTEXT_H
#define	CEC_CONTEXT_H

#include "matrix.h"
#include "model.h"
#include "alloc.h"

struct cec_matrix_array {
    int len;
    struct cec_matrix * mats[];
};

typedef struct cec_matrix_array cec_mat_arr;

struct cec_temp_data {
    struct cec_matrix * t_mean_matrix;
    struct cec_matrix * t_matrix_nn;
    struct cec_matrix * n_covariance_matrix;
    struct cec_matrix_array  * t_covariance_matrices;
};

struct cec_input {
    const struct cec_matrix * points;
    const struct cec_matrix * centers;
    const vec_i * initial_assignment;
    const cec_model * const* models;
    int max_iterations;
    int min_card;
};

struct cec_output {
    cec_mat * centers;
    vec_i * clustering_vector;
    vec_i * clusters_number;
    vec_d * energy;
    int iterations;
    struct cec_matrix_array * covriances;
};

typedef struct cec_output cec_out;
typedef struct cec_input cec_in;
typedef struct cec_temp_data cec_tmp_data;

struct cec_matrix_array * create_cec_matrix_array(int len, int m, int n);

cec_out * create_cec_output(int m, int k, int n, int max_iterations);

cec_in * create_cec_input(const cec_mat * points, const cec_mat * centers, const vec_i * initial_assignment,
                                    const cec_model ** models, int max_iterations, int min_card);
cec_tmp_data * create_temp_data(int k, int n);

void cec_copy_results_content(const cec_out *from, cec_out *to);

double cec_final_energy(cec_out * c_res);
int cec_final_centers_number(cec_out *c_res);

#endif	/* CEC_CONTEXT_H */
