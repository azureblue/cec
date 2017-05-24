#ifndef CEC_CONTEXT_H
#define	CEC_CONTEXT_H

#include "matrix.h"
#include "energy.h"
#include "model.h"
#include "alloc.h"

struct cec_matrix_array {
    int l;
    struct cec_matrix * mats[];
};

struct cec_temp_data
{
    struct cec_matrix * t_mean_matrix;
    struct cec_matrix * t_matrix_nn;
    struct cec_matrix * n_covariance_matrix;
    struct cec_matrix_array  * t_covariance_matrices;
};

struct cec_input
{
    const struct cec_matrix * points;
    const struct cec_matrix * centers;
    struct cec_model ** models;
    int max_iterations;
    int min_card;
};

struct cec_results
{
    struct cec_matrix * centers;
    int * clustering_vector;
    int * clusters_number;
    int iterations;
    double * energy;
    struct cec_matrix_array * covriances;
    int error;
};

struct cec_context
{
    struct cec_input * input;
    struct cec_results * results;
    struct cec_temp_data * temp_data;
};

typedef struct cec_results cec_res;
typedef struct cec_input cec_in;
typedef struct cec_context cec_ctx;

struct cec_matrix_array * create_cec_matrix_array(int l, int m, int n);

struct cec_context * create_cec_context(struct cec_input *in, struct cec_results *res);

struct cec_results * create_cec_result(int m, int k, int n, int max_iterations);

struct cec_input * create_cec_input(const struct cec_matrix * points, const struct cec_matrix * centers,
                                    struct cec_model ** models, int max_iterations, int min_card);

void cec_copy_results_content(cec_res *from, cec_res *to, int m);

double cec_final_energy(cec_res * c_res);

#endif	/* CEC_CONTEXT_H */
