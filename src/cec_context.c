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

struct cec_results * create_cec_result(int m, int k, int n, int max_iterations) {
    struct cec_results *results = alloc(struct cec_results);
    results->error = NO_ERROR;
    results->clustering_vector = alloc_n(int, m);
    results->centers = cec_matrix_create(k, n);
    results->clusters_number = alloc_n(int, max_iterations + 1);
    results->energy = alloc_n(double, max_iterations + 1);
    results->covriances = create_cec_matrix_array(k, n, n);
    return results;
}

struct cec_input * create_cec_input(const struct cec_matrix * points, const struct cec_matrix * centers,
                                    struct cec_model ** models, int max_iterations, int min_card) {
    struct cec_input *input = alloc(struct cec_input);
    input->points = points;
    input->centers = centers;
    input->models = models;
    input->max_iterations = max_iterations;
    input->min_card = min_card;
    return input;
}

struct cec_context * create_cec_context(struct cec_input *in, struct cec_results *res)
{
    int k = in->centers->m;
    int n = in->points->n;
    struct cec_context * context = alloc(struct cec_context);
    context->input = in;
    context->results = res;
    context->temp_data = create_temp_data(k, n);
    
    return context;
}

double cec_final_energy(cec_res *c_res) {
    return c_res->energy[c_res->iterations];
}

void cec_copy_results_content(cec_res *from, cec_res *to, int m) {
    to->iterations = from->iterations;
    int output_size = to->iterations + 1;
    for (int i = 0; i < output_size; i++) {
        to->energy[i] = from->energy[i];
        to->clusters_number[i] = from->clusters_number[i];
    }
    for (int i = 0; i < m; i++)
        to->clustering_vector[i] = from->clustering_vector[i];
    cec_matrix_copy_data(from->centers, to->centers);
    for (int i = 0; i < from->covriances->l; i++)
        cec_matrix_copy_data(from->covriances->mats[i], to->covriances->mats[i]);
}
