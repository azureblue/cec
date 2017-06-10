#include "energy.h"
#include "../cov_utils.h"
#include "../alloc.h"
#include "../model.h"

static const int EIGENVALUES_WORKSPACE_BLOCK = 128;

struct context_fe
{
    double * evals_given;
    double * evals_temp;
    double given_evals_product;
    struct cec_matrix * workspace;
    struct cec_matrix * temp_matrix;
};

static double cross_entropy(struct cross_entropy_context * context, const struct cec_matrix * cov) {
    struct context_fe * cfe = (struct context_fe *) context->specific_context;
    struct cec_matrix * temp_matrix = cfe->temp_matrix;
    double * given_evals = cfe -> evals_given;
    double * evals = cfe -> evals_temp;
    int n = cov->n;
    int error = cec_cov_eigenvalues(cov, temp_matrix, cfe->workspace, evals);
    if (error) {
        context->last_error = error;
        return NAN;
    }
    double e_sum = 0;
    for (int i = 0; i < n; i++)
        e_sum += evals[i] / given_evals[i];
    return (n / 2.0) * log(2.0 * M_PI) + (1.0 / 2.0) * e_sum + (1.0 / 2.0)
                                                               * log(cfe->given_evals_product);
}

struct cec_model * cec_create_ce_ctx_eigenvalues(int n, const double * given_evals) {
    struct cec_model *model = alloc(struct cec_model);
    model->cross_entropy_context = alloc(cross_entropy_ctx);
    struct context_fe * c_fe =  alloc(struct context_fe);
    c_fe->evals_given = alloc_n(double, n);
    array_copy(given_evals, c_fe->evals_given, n);
    c_fe->evals_temp = alloc_n(double, n);
    c_fe->workspace = cec_matrix_create(EIGENVALUES_WORKSPACE_BLOCK, 1);
    c_fe->temp_matrix = cec_matrix_create(n, n);
    double product = 1.0;
    for (int i = 0; i < n; i++)
        product *= given_evals[i];
    c_fe->given_evals_product = product;
    model->cross_entropy_context->specific_context = c_fe;
    model->cross_entropy = cross_entropy;
    return model;
}
