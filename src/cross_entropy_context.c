#include <stdlib.h>

#include "alloc.h"
#include "cross_entropy_context.h"
#include "matrix.h"

static const int EIGENVALUES_WORKSPACE_BLOCK = 128;

void destroy_cross_entropy_context(struct cross_entropy_context * context)
{
    if (context == NULL)
	return;

    cec_matrix_destroy(context->temp_matrix);
    switch (context->family)
    {
	case FIXED_R:
	    m_free(context->custom_context);
	    break;
	case GIVEN_COVARIANCE:
	{
	    struct context_gc * c_gc = (struct context_gc *) context->custom_context;
	    if (context->custom_context == NULL)
		break;

	    cec_matrix_destroy(c_gc->given_cov);
	    cec_matrix_destroy(c_gc->i_given_cov);
	    m_free(context->custom_context);
	    break;
	}
	case FIXEDEIGENVALUES:
	{
	    if (context->custom_context == NULL)
		break;

	    struct context_fe * c_fe = (struct context_fe *) context->custom_context;
	    m_free(c_fe->evals);
	    m_free(c_fe->given_evals);
	    cec_matrix_destroy(c_fe->workspace);
	    m_free(context->custom_context);
	    break;

	}
	case ALL:
	case SPHERICAL:
	case DIAGONAL:
	    break;
    }

    m_free(context);
}

struct cross_entropy_context * create_cross_entropy_context(
	enum density_family family, int n)
{
    struct cross_entropy_context * context = m_alloc(
	    sizeof (struct cross_entropy_context));
    if (context == NULL)
	return NULL;

    context->n = n;
    context->custom_context = NULL;
    context->temp_matrix = cec_matrix_create(n, n);
    context->family = family;

    if (context->temp_matrix == NULL)
    {
	destroy_cross_entropy_context(context);
	return NULL;
    }

    switch (family)
    {
	case FIXED_R:
	{
	    struct context_r * c_r = m_alloc(sizeof (struct context_r));
	    if (c_r == NULL)
	    {
		destroy_cross_entropy_context(context);
		return NULL;
	    }
	    context->custom_context = c_r;
	    break;
	}
	case GIVEN_COVARIANCE:
	{
	    struct context_gc * context_gc = m_alloc(sizeof (struct context_gc));
	    if (context_gc == NULL)
	    {
		destroy_cross_entropy_context(context);
		return NULL;
	    }
	    context->custom_context = context_gc;

	    context_gc->given_cov = cec_matrix_create(n, n);
	    context_gc->i_given_cov = cec_matrix_create(n, n);
	    if (context_gc->given_cov == NULL || context_gc->i_given_cov == NULL)
	    {
		destroy_cross_entropy_context(context);
		return NULL;
	    }
	    break;
	}
	case FIXEDEIGENVALUES:
	{
	    struct context_fe * context_fe = m_alloc(sizeof (struct context_fe));
	    if (context_fe == NULL)
	    {
		destroy_cross_entropy_context(context);
		return NULL;
	    }
	    context->custom_context = context_fe;
	    
	    context_fe->evals	    = m_alloc(sizeof (double) * n);
	    context_fe->given_evals = m_alloc(sizeof (double) * n);
	    context_fe->workspace   = cec_matrix_create((EIGENVALUES_WORKSPACE_BLOCK + 2) * n, 1);

	    if (context_fe->evals == NULL || context_fe->given_evals == NULL || context_fe->workspace == NULL)
	    {
		destroy_cross_entropy_context(context);
		return NULL;
	    }
	    break;
	}

	case ALL:
	case SPHERICAL:
	case DIAGONAL:
	    break;
    }

    return context;
}
