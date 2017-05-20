#include "alloc.h"
#include "cec_context.h"

struct cec_matrix_array * create_cec_matrix_array(int l, int m, int n)
{
    struct cec_matrix_array * matrix_array = alloc_fam(struct cec_matrix_array, struct cec_matrix, l);
    
    matrix_array->l = l;
        
    for (int i = 0; i < l; i++)
        matrix_array->mats[i] = cec_matrix_create(m, n);
    
    return matrix_array;
}

static struct cec_temp_data * create_temp_data(int k, int n)
{
    struct cec_temp_data * data = alloc(struct cec_temp_data);

    data->n_covariance_matrix = cec_matrix_create(n, n);
    data->t_matrix_nn = cec_matrix_create(n, n);
    data->t_mean_matrix = cec_matrix_create(k, n);
    data->t_covariance_matrices = create_cec_matrix_array(k, n, n);

    return data;
}

struct cec_context * create_cec_context(const struct cec_matrix * points, const struct cec_matrix * centers,
        struct cec_model ** models, int max_iterations, int min_card)
{
    int m = points->m;
    int k = centers->m;
    int n = points->n;
    
    struct cec_context * context = alloc(struct cec_context);
    context->input = alloc(struct cec_input);
    context->results = alloc(struct cec_results);
    
    context->input->points = points;
    context->input->centers = centers;
    context->input->models = models;
    context->input->max_iterations = max_iterations;
    context->input->min_card = min_card;    
    context->results->error = NO_ERROR;
    context->results->clustering_vector = alloc_n(int, m);
    context->results->centers = cec_matrix_create(k, n);
    context->results->clusters_number = alloc_n(int, max_iterations + 1);
    context->results->energy = alloc_n(double, max_iterations + 1);
    context->results->covriances = create_cec_matrix_array(k, n, n);
    context->temp_data = create_temp_data(k, n);
    
    return context;
}
