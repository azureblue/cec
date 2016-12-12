#include "alloc.h"
#include "alloc_check.h"
#include "cross_entropy_context.h"
#include "matrix.h"
#include "model.h"

static const int EIGENVALUES_WORKSPACE_BLOCK = 128;
 
struct cross_entropy_context * create_cross_entropy_context_all(int n)
{
    checked_alloc(3);
    
    struct cross_entropy_context * context = check_alloc(struct cross_entropy_context);
    struct context_all * c_all = check_alloc(struct context_all);
    
    c_all->temp_matrix = check_ptr(cec_matrix_create(n, n));
    
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
    checked_alloc(2);
    
    struct cross_entropy_context * context = check_alloc(struct cross_entropy_context);
    struct context_r * c_r =  check_alloc(struct context_r);
    
    c_r->r = r;
    
    context->specific_context = c_r;
    
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_covariance(int n, const struct cec_matrix * cov, const struct cec_matrix * cov_i)
{
    checked_alloc(5);
    
    struct cross_entropy_context * context = check_alloc(struct cross_entropy_context);
    struct context_gc * c_gc =  check_alloc(struct context_gc);
    
    c_gc->given_cov = check_ptr(cec_matrix_create_copy(cov));
    c_gc->i_given_cov = check_ptr(cec_matrix_create_copy(cov_i));
    c_gc->temp_matrix = check_ptr(cec_matrix_create(n, n));
    
    context->specific_context = c_gc;
    
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_eigenvalues(int n, const double * given_evals)
{
    checked_alloc(6);
    
    struct cross_entropy_context * context = check_alloc(struct cross_entropy_context);     
    struct context_fe * c_fe =  check_alloc(struct context_fe);
    
    c_fe->evals_given = check_alloc_n(double, n);
    array_copy(given_evals, c_fe->evals_given, n);
    
    c_fe->evals_temp = check_alloc_n(double, n);
    c_fe->workspace = check_ptr(cec_matrix_create(EIGENVALUES_WORKSPACE_BLOCK, 1));
    c_fe->temp_matrix = check_ptr(cec_matrix_create(n, n));
    
    double product = 1.0;
    for (int i = 0; i < n; i++)
        product *= given_evals[i];
    
    c_fe->given_evals_product = product;
    
    context->specific_context = c_fe;    
    
    return context;
}

void destroy_cross_entropy_context_eigenvalues(struct cross_entropy_context * context)
{
    if (!context) 
        return;
    
    struct context_fe * specific_context = (struct context_fe *) context->specific_context;
    
    m_free(specific_context->evals_given);
    m_free(specific_context->evals_temp);
    m_free(specific_context->temp_matrix);
    m_free(specific_context->workspace);
    m_free(specific_context);
    m_free(context);
}

void destroy_cross_entropy_context_covariance(struct cross_entropy_context * context)
{
    if (!context) 
        return;
    
    struct context_gc * specific_context = (struct context_gc *) context->specific_context;
    
    m_free(specific_context->given_cov);
    m_free(specific_context->i_given_cov);
    m_free(specific_context->temp_matrix);
    m_free(specific_context);
    m_free(context);
}

void destroy_cross_entropy_context_fixedr(struct cross_entropy_context * context)
{
    if (!context) 
        return;
    
    struct context_r * specific_context = (struct context_r *) context->specific_context;
    
    m_free(specific_context);
    m_free(context);
}

void destroy_cross_entropy_context_all(struct cross_entropy_context * context)
{
    if (!context) 
        return;
    
    struct context_all * specific_context = (struct context_all *) context->specific_context;
    
    m_free(specific_context->temp_matrix);
    m_free(specific_context);
    m_free(context);
}

void destroy_cross_entropy_context_diagonal(struct cross_entropy_context * context)
{   
    m_free(context);
}

void destroy_cross_entropy_context_spherical(struct cross_entropy_context * context)
{  
    m_free(context);
}
