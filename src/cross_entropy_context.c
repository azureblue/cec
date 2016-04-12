#include <stdlib.h>

#include "alloc.h"
#include "cross_entropy_context.h"
#include "matrix.h"
#include "model.h"
#include "checked_allocator.h"

static const int EIGENVALUES_WORKSPACE_BLOCK = 128;
 
struct cross_entropy_context * create_cross_entropy_context_all(int n)
{
    checked_allocation;
    
    struct cross_entropy_context * context = check(m_alloc(sizeof (struct cross_entropy_context)));
    struct context_all * c_all = check(m_alloc(sizeof (struct context_all)));
    c_all->temp_matrix = check(cec_matrix_create(n, n));
    
    context->specific_context = c_all;
    
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_spherical(int n)
{
  return m_alloc(sizeof (struct cross_entropy_context));
}

struct cross_entropy_context * create_cross_entropy_context_diagonal(int n)
{
    return m_alloc(sizeof (struct cross_entropy_context));
}

struct cross_entropy_context * create_cross_entropy_context_fixedr(int n, double r)
{
    checked_allocation;
    
    struct cross_entropy_context * context = check(m_alloc(sizeof (struct cross_entropy_context)));
    struct context_r * c_r =  check(m_alloc(sizeof (struct context_r)));
    c_r->r = r;
    
    context->specific_context = c_r;
    
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_covariance(int n, const struct cec_matrix * cov, const struct cec_matrix * cov_i)
{
    checked_allocation;
    
    struct cross_entropy_context * context = check(m_alloc(sizeof (struct cross_entropy_context)));
    struct context_gc * c_gc =  check(m_alloc(sizeof (struct context_gc)));
    c_gc->given_cov = check(cec_matrix_create_copy(cov));
    c_gc->i_given_cov = check(cec_matrix_create_copy(cov_i));
    c_gc->temp_matrix = check(cec_matrix_create(n, n));
    
    context->specific_context = c_gc;
    
    return context;
}

struct cross_entropy_context * create_cross_entropy_context_eigenvalues(int n, const double * given_evals)
{
    checked_allocation;
    
    struct cross_entropy_context * context = check(m_alloc(sizeof (struct cross_entropy_context)));
     
    struct context_fe * c_fe =  check(m_alloc(sizeof (struct context_fe)));
    
    c_fe->evals_given = check(m_alloc(sizeof (double) * n));
    array_copy(given_evals, c_fe->evals_given, n);
    
    c_fe->evals_temp = check(m_alloc(sizeof (double) * n));
    c_fe->workspace = check(cec_matrix_create(EIGENVALUES_WORKSPACE_BLOCK, 1));
    c_fe->temp_matrix = check(cec_matrix_create(n, n));
    
    
    double p = 1;
    for (int i = 0; i < n; i++)
        p *= given_evals[i];
    
    c_fe->given_evals_product = p;
    
    context->specific_context = c_fe;    
    
    return context;
}

void destroy_cross_entropy_context_eigenvalues(struct cross_entropy_context * context)
{
    if (!context) return;
    
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
    if (!context) return;
    
    struct context_gc * specific_context = (struct context_gc *) context->specific_context;
    m_free(specific_context->given_cov);
    m_free(specific_context->i_given_cov);
    m_free(specific_context->temp_matrix);
    m_free(specific_context);
    m_free(context);
}

void destroy_cross_entropy_context_fixedr(struct cross_entropy_context * context)
{
    if (!context) return;
    
    struct context_r * specific_context = (struct context_r *) context->specific_context;
    m_free(specific_context);
    m_free(context);
}

void destroy_cross_entropy_context_all(struct cross_entropy_context * context)
{
    if (!context) return;
    
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
