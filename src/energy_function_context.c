#include <stdlib.h>

#include "alloc.h"
#include "energy_function_context.h"
#include "matrix.h"

void destroy_energy_function_context(struct energy_function_context * context)
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
		break;
	}
	case ALL:
	case SPHERICAL:
	case DIAGONAL:
		break;
	}

	m_free(context);
}

struct energy_function_context * create_energy_function_context(
		enum density_family family, int n)
{
	struct energy_function_context * context = m_alloc(
			sizeof(struct energy_function_context));
	if (context == NULL)
		return NULL;

	context->n = n;
	context->custom_context = NULL;
	context->temp_matrix = cec_matrix_create(n, n);
	context->family = family;

	if (context->temp_matrix == NULL)
	{
		destroy_energy_function_context(context);
		return NULL;
	}

	switch (family)
	{
	case FIXED_R:
	{
		struct context_r * c_r = m_alloc(sizeof(struct context_r));
		if (c_r == NULL)
		{
			destroy_energy_function_context(context);
			return NULL;
		}
		context->custom_context = c_r;
		break;
	}
	case GIVEN_COVARIANCE:
	{
		struct context_gc * context_gc = m_alloc(sizeof(struct context_gc));
		if (context_gc == NULL)
		{
			destroy_energy_function_context(context);
			return NULL;
		}

		context_gc->given_cov = cec_matrix_create(n, n);
		context_gc->i_given_cov = cec_matrix_create(n, n);
		if (context_gc->given_cov == NULL || context_gc->i_given_cov == NULL)
		{
			destroy_energy_function_context(context);
			return NULL;
		}
		context->custom_context = context_gc;
		break;
	}
	case ALL:
	case SPHERICAL:
	case DIAGONAL:
		break;
	}

	return context;
}
