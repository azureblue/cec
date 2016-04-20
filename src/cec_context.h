#ifndef CEC_CONTEXT_H
#define	CEC_CONTEXT_H

#include "matrix.h"
#include "energy.h"
#include "model.h"

struct cec_temp_data
{
    struct cec_matrix * t_mean_matrix;
    struct cec_matrix * t_matrix_nn;
    struct cec_matrix * n_covariance_matrix;
    struct cec_matrix ** t_covariance_matrices;
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
    struct cec_matrix ** covriances;
    int error;
};

struct cec_context
{
    /*
     * Input parameters.
     */
    struct cec_input * input;

    /*
     * CEC result (output parameters).
     * Memory must be allocated before performing the algorithm.
     */
    struct cec_results * results;

    /*
     * Temporary data. 
     * Memory must be allocated before performing the algorithm.
     */
    struct cec_temp_data * temp_data;
};

struct cec_matrix ** create_cec_matrix_array(int m, int n, int l);

struct cec_context * create_cec_context(
        const struct cec_matrix * points,
        const struct cec_matrix * centers,
        struct cec_model **,
        int max_iterations,
        int min_card
        );

void destroy_cec_matrix_array(struct cec_matrix ** matrix_array, int l);

void destroy_cec_context_results(struct cec_context * context);

#endif	/* CEC_CONTEXT_H */

