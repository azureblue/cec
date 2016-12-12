#include "alloc.h"
#include "alloc_check.h"
#include "cec_context.h"

void destroy_cec_matrix_array(memptr_t matrix_array_ptr)
{
    struct cec_matrix_array * matrix_array = matrix_array_ptr;
    
    if (!matrix_array)
        return;

    int l = matrix_array->l;
    
    for (int i = 0; i < l; i++)
        cec_matrix_destroy(matrix_array->mats[i]);

    m_free(matrix_array);
}

static void destroy_cec_temp_data(memptr_t data_ptr)
{
    struct cec_temp_data * data = data_ptr;
    
    if (!data)
        return;

    cec_matrix_destroy(data->n_covariance_matrix);
    cec_matrix_destroy(data->t_matrix_nn);
    cec_matrix_destroy(data->t_mean_matrix);

    destroy_cec_matrix_array(data->t_covariance_matrices);

    m_free(data);
}

void destroy_cec_context_results(struct cec_context * context)
{
    if (!context)
        return;

    destroy_cec_matrix_array(context->results->covriances);
    destroy_cec_temp_data(context->temp_data);

    m_free(context->results->clustering_vector);
    m_free(context->results->clusters_number);
    m_free(context->results->energy);
    m_free(context->results->centers);
    m_free(context->results);
    m_free(context->input);
    m_free(context);
}

struct cec_matrix_array * create_cec_matrix_array(int l, int m, int n)
{
    struct cec_matrix_array * matrix_array = m_alloc(
            sizeof (struct cec_matrix_array) + sizeof (struct cec_matrix *) * l);
    
    if (!matrix_array)
        return NULL;
    
    matrix_array->l = l;
    
    for (int i = 0; i < l; i++)
        matrix_array->mats[i] = NULL;
    
    for (int i = 0; i < l; i++)
        if (!(matrix_array->mats[i] = cec_matrix_create(m, n)))
        {
            destroy_cec_matrix_array(matrix_array);
            return NULL;
        }

    return matrix_array;
}

static struct cec_temp_data * create_temp_data(int k, int n)
{
    checked_alloc(5);
    
    struct cec_temp_data * data = check_alloc(struct cec_temp_data);

    data->n_covariance_matrix = check_ptr(cec_matrix_create(n, n));
    data->t_matrix_nn = check_ptr(cec_matrix_create(n, n));
    data->t_mean_matrix = check_ptr(cec_matrix_create(k, n));
    data->t_covariance_matrices = check_ptr_dstr(create_cec_matrix_array(k, n, n), destroy_cec_matrix_array);

    return data;
}

struct cec_context *
create_cec_context(const struct cec_matrix * points, const struct cec_matrix * centers,
        struct cec_model ** models, int max_iterations, int min_card)
{
    checked_alloc(9);
    
    struct cec_context * context = check_alloc(struct cec_context);
    context->input = check_alloc(struct cec_input);
    context->results = check_alloc(struct cec_results);
    
    int m = points->m;
    int k = centers->m;
    int n = points->n;

    context->input->points = points;
    context->input->centers = centers;
    context->input->models = models;
    context->input->max_iterations = max_iterations;
    context->input->min_card = min_card;
    
    context->results->error = NO_ERROR;
    context->results->clustering_vector = check_alloc_n(int, m);
    context->results->centers = check_ptr(cec_matrix_create(k, n));
    context->results->clusters_number = check_alloc_n(int, max_iterations + 1);
    context->results->energy = check_alloc_n(double, max_iterations + 1);
    context->results->covriances = check_ptr_dstr(create_cec_matrix_array(k, n, n), destroy_cec_matrix_array);
    
    context->temp_data = check_ptr_dstr(create_temp_data(k, n), destroy_cec_temp_data);
    
    return context;
}
