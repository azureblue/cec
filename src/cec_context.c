#include "alloc.h"
#include "cec_context.h"

cec_mat_arr * create_cec_matrix_array(int len, int m, int n)
{
    cec_mat_arr * matrix_array = alloc_fam(cec_mat_arr, cec_mat*, len);
    matrix_array->len = len;
    for (int i = 0; i < len; i++)
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

cec_out * create_cec_output(int m, int k, int n, int max_iterations) {
    cec_out *output = alloc(cec_out);
    output->clustering_vector = vec_i_create(m);
    output->clusters_number = vec_i_create(max_iterations + 1);
    output->energy = vec_d_create(max_iterations + 1);
    output->centers = cec_matrix_create(k, n);
    output->covriances = create_cec_matrix_array(k, n, n);
    output->initial_k = k;
    return output;
}

cec_in * create_cec_input(const cec_mat * points, cec_mat * centers,
                                    struct cec_model ** models, int max_iterations, int min_card) {
    cec_in *input = alloc(cec_in);
    input->points = points;
    input->centers = centers;
    input->models = models;
    input->max_iterations = max_iterations;
    input->min_card = min_card;
    return input;
}

cec_ctx * create_cec_context(cec_in *in, cec_out *out)
{    
    int k = in->centers->m;
    int n = in->points->n;
    cec_ctx * ctx = alloc(cec_ctx);
    ctx->input = in;
    ctx->results = out;
    ctx->temp_data = create_temp_data(k, n);
    
    return ctx;
}

double cec_final_energy(cec_out *c_res) {
    return c_res->energy->ar[c_res->iterations];
}

void cec_copy_results_content(cec_out *from, cec_out *to) {
    to->initial_k = from->initial_k;
    to->iterations = from->iterations;
    vec_i_copy(from->clustering_vector, to->clustering_vector);
    vec_i_copy(from->clusters_number, to->clusters_number);
    vec_d_copy(from->energy, to->energy);
    cec_matrix_copy_data(from->centers, to->centers);
    for (int i = 0; i < from->covriances->len; i++)
        cec_matrix_copy_data(from->covriances->mats[i], to->covriances->mats[i]);
}
