#include "cec_context.h"
#include "alloc.h"
#include "checked_allocator.h"

void destroy_cec_matrix_array(struct cec_matrix ** matrix_array, int l)
{
    if (!matrix_array)
        return;

    for (int i = 0; i < l; i++)
        cec_matrix_destroy(matrix_array[i]);

    m_free(matrix_array);
}

static void destroy_cec_temp_data(struct cec_temp_data * data, int k)
{
    if (!data)
        return;

    cec_matrix_destroy(data->n_covariance_matrix);
    cec_matrix_destroy(data->t_matrix_nn);
    cec_matrix_destroy(data->t_mean_matrix);

    destroy_cec_matrix_array(data->t_covariance_matrices, k);

    m_free(data);
}

void destroy_cec_context_results(struct cec_context * context)
{
    if (!context)
        return;

    int k = context->input->centers->m;

    destroy_cec_matrix_array(context->results->covriances, k);
    destroy_cec_temp_data(context->temp_data, k);

    m_free(context->results->clustering_vector);
    m_free(context->results->clusters_number);
    m_free(context->results->energy);
    m_free(context->results->centers);
    m_free(context->results);
    m_free(context->input);
    m_free(context);
}

struct cec_matrix ** create_cec_matrix_array(int l, int m, int n)
{
    struct cec_matrix ** matrix_array = m_alloc(
            sizeof (struct cec_matrix *) * l);
    
    if (!matrix_array)
        return NULL;

    int m_alloc_error_flag = 0;

    for (int i = 0; i < l; i++)
    {
        matrix_array[i] = cec_matrix_create(m, n);
        if (!matrix_array[i])
            m_alloc_error_flag = 1;
    }

    if (m_alloc_error_flag)
    {
        destroy_cec_matrix_array(matrix_array, l);
        return NULL;
    }

    return matrix_array;
}

static struct cec_temp_data * create_temp_data(int k, int n)
{
    checked_allocation;
    struct cec_temp_data * data = check(m_alloc(sizeof (struct cec_temp_data)));

    data->n_covariance_matrix = check(cec_matrix_create(n, n));
    data->t_matrix_nn = check(cec_matrix_create(n, n));
    data->t_mean_matrix = check(cec_matrix_create(k, n));
    if (!(data->t_covariance_matrices = create_cec_matrix_array(k, n, n)))
    {
        free_checked_pointers();
        return NULL;
    }

    return data;
}

struct cec_context *
create_cec_context(const struct cec_matrix * points, const struct cec_matrix * centers,
        struct cec_model ** models, int max_iterations, int min_card)
{
    checked_allocation;
    struct cec_context * context = check(m_alloc(sizeof (struct cec_context)));
    context->input = check(m_alloc(sizeof (struct cec_input)));
    context->results = check(m_alloc(sizeof (struct cec_results)));
    
    int m = points->m;
    int k = centers->m;
    int n = points->n;

    context->input->points = points;
    context->input->centers = centers;
    context->input->models = models;
    context->input->max_iterations = max_iterations;
    context->input->min_card = min_card;
    
    context->results->error = NO_ERROR;
    context->results->clustering_vector = check(m_alloc(sizeof (int) * (m)));
    context->results->centers = check(cec_matrix_create(k, n));
    context->results->clusters_number = check(m_alloc(sizeof (int) * (max_iterations + 1)));
    context->results->energy = check(m_alloc(sizeof (double) * (max_iterations + 1)));

    if (!(context->results->covriances = create_cec_matrix_array(k, n, n)))
    {
        free_checked_pointers();
        return NULL;
    }
    
    if (!(context->temp_data = create_temp_data(k, n)))
    {
        destroy_cec_matrix_array(context->results->covriances, k);
        free_checked_pointers();
        return NULL;
    }
    
    return context;
}
