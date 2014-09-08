#include "cec_context.h"
#include "alloc.h"

void destroy_cec_matrix_array(struct cec_matrix ** matrix_array, int l)
{
    if (matrix_array == NULL)
	return;

    for (int i = 0; i < l; i++)
	cec_matrix_destroy(matrix_array[i]);

    m_free(matrix_array);
}

static void destroy_cec_temp_data(struct cec_temp_data * data, int k)
{
    if (data == NULL)
	return;

    cec_matrix_destroy(data->n_covariance_matrix);
    cec_matrix_destroy(data->t_matrix_nn);
    cec_matrix_destroy(data->t_mean_matrix);

    destroy_cec_matrix_array(data->t_covariance_matrices, k);

    m_free(data);
}

void destroy_cec_context_results(struct cec_context * context)
{
    if (context == NULL)
	return;

    int k = context->centers->m;

    destroy_cec_matrix_array(context->covriances, k);
    destroy_cec_temp_data(context->temp_data, k);

    m_free(context->clustering_vector);
    m_free(context->clusters_number);
    m_free(context->energy);
    m_free(context);
}

struct cec_matrix ** create_cec_matrix_array(int l, int m, int n)
{
    struct cec_matrix ** matrix_array = m_alloc(
	    sizeof (struct cec_matrix *) * l);
    if (matrix_array == NULL)
	return NULL;

    int m_alloc_error_flag = 0;

    for (int i = 0; i < l; i++)
    {
	matrix_array[i] = cec_matrix_create(m, n);
	if (matrix_array[i] == NULL)
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
    struct cec_temp_data * data = m_alloc(sizeof (struct cec_temp_data));
    if (data == NULL)
	return NULL;

    data->n_covariance_matrix = cec_matrix_create(n, n);
    data->t_matrix_nn = cec_matrix_create(n, n);
    data->t_mean_matrix = cec_matrix_create(k, n);
    data->t_covariance_matrices = create_cec_matrix_array(k, n, n);

    if (data->n_covariance_matrix == NULL || data->t_matrix_nn == NULL
	    || data->t_mean_matrix == NULL
	    || data->t_covariance_matrices == NULL)
    {
	destroy_cec_temp_data(data, k);
	return NULL;
    }

    return data;
}

struct cec_context *
create_cec_context(struct cec_matrix * points, struct cec_matrix * centers,
	struct energy_function_context ** energy_function_contexts,
	energy_function * energy_functions, int max_iterations, int min_card)
{
    struct cec_context * context = m_alloc(sizeof (struct cec_context));

    if (context == NULL)
	return NULL;

    int m = points->m;
    int k = centers->m;
    int n = points->n;

    context->points = points;
    context->centers = centers;
    context->energy_function_contexts = energy_function_contexts;
    context->energy_functions = energy_functions;
    context->max_iterations = max_iterations;
    context->min_card = min_card;
    context->error = NO_ERROR;

    context->clustering_vector = m_alloc(sizeof (int) * (m));
    context->clusters_number = m_alloc(sizeof (int) * (max_iterations + 1));
    context->energy = m_alloc(sizeof (double) * (max_iterations + 1));

    context->covriances = create_cec_matrix_array(k, n, n);
    context->temp_data = create_temp_data(k, n);

    if (context->clustering_vector == NULL
	    || (context->energy == NULL && max_iterations >= 0)
	    || (context->clusters_number == NULL  && max_iterations >= 0)
	    || context->covriances == NULL
	    || context->temp_data == NULL)
    {
	destroy_cec_context_results(context);
	return NULL;
    }

    return context;
}
