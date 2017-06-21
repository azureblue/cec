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
    struct cec_model ** models;
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
    int initial_k;
};

struct cec_context
{
    struct cec_input * input;
    struct cec_output * results;
    struct cec_temp_data * temp_data;
};

typedef struct cec_output cec_out;
typedef struct cec_input cec_in;
typedef struct cec_context cec_ctx;

struct cec_matrix_array * create_cec_matrix_array(int len, int m, int n);

cec_ctx * create_cec_context(cec_in *in, cec_out *out);

cec_out * create_cec_output(int m, int k, int n, int max_iterations);

cec_in * create_cec_input(const cec_mat * points, const cec_mat * centers,
                                    struct cec_model ** models, int max_iterations, int min_card);

void cec_copy_results_content(cec_out *from, cec_out *to);

double cec_final_energy(cec_out * c_res);

#endif	/* CEC_CONTEXT_H */
