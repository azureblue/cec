#include "alloc.h"
#include "cross_entropy_context.h"
#include "matrix.h"

static const int EIGENVALUES_WORKSPACE_BLOCK = 128;
 
struct cross_entropy_context * create_cross_entropy_context_all(int n)
{
    struct cross_entropy_context * context = alloc(struct cross_entropy_context);
    struct context_all * c_all = alloc(struct context_all);    
    c_all->temp_matrix = cec_matrix_create(n, n);    
    context->specific_context = c_all;    
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_spherical()
{
    return alloc(struct cross_entropy_context);
}

struct cross_entropy_context * create_cross_entropy_context_diagonal()
{
    return alloc(struct cross_entropy_context);
}

struct cross_entropy_context * create_cross_entropy_context_fixedr(double r)
{
    struct cross_entropy_context * context = alloc(struct cross_entropy_context);
    struct context_r * c_r =  alloc(struct context_r);    
    c_r->r = r;    
    context->specific_context = c_r;    
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_covariance(int n, const struct cec_matrix * cov, const struct cec_matrix * cov_i)
{
    struct cross_entropy_context * context = alloc(struct cross_entropy_context);
    struct context_gc * c_gc =  alloc(struct context_gc);    
    c_gc->given_cov = cec_matrix_create_copy(cov);
    c_gc->i_given_cov = cec_matrix_create_copy(cov_i);
    c_gc->temp_matrix = cec_matrix_create(n, n);
    context->specific_context = c_gc;
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_eigenvalues(int n, const double * given_evals)
{
    struct cross_entropy_context * context = alloc(struct cross_entropy_context);     
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
    context->specific_context = c_fe;    
    return context;
}
