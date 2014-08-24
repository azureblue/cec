#include "energy_function.h"

energy_function energy_function_for(enum density_family family)
{
	switch (family)
	{
	case ALL:
		return h_all;
	case SPHERICAL:
		return h_spherical;
	case DIAGONAL:
		return h_diagonal;
	case FIXED_R:
		return h_fixed_r;
	case GIVEN_COVARIANCE:
		return h_given_covariance;
	}

	return NULL;
}

double h_given_covariance(const struct energy_function_context * context,
		const struct cec_matrix * cov)
{
	struct context_gc * cgc = (struct context_gc *) context->custom_context;
	struct cec_matrix * temp_matrix = context->temp_matrix;
	int n = context->n;

	cec_matrix_multiply_sq(cgc->i_given_cov, cov, temp_matrix);
	double trace = cec_matrix_trace(temp_matrix);

	if (isnan(trace))
	{
		*context->last_error = POSITIVE_DEFINITE_ERROR;
		return NAN;
	}

	double det = cec_matrix_det_positive_definite(cgc->given_cov, temp_matrix);
	if (isnan(det))
	{
		*context->last_error = POSITIVE_DEFINITE_ERROR;
		return NAN;
	}

	return (n / 2.0) * log(2.0 * M_PI) + (1.0 / 2.0) * trace
			+ (1.0 / 2.0) * log(det);
}

double h_spherical(const struct energy_function_context * context,
		const struct cec_matrix * cov)
{
	int n = context->n;

	double trace = cec_matrix_trace_assert_positive(cov);
	if (isnan(trace))
	{
		*context->last_error = POSITIVE_DEFINITE_ERROR;
		return NAN;
	}
	return (n / 2.0) * log(2.0 * M_PI * M_E / n) + (n / 2.0) * log(trace);
}

double h_fixed_r(const struct energy_function_context * context,
		const struct cec_matrix * cov)
{
	struct context_r * cr = (struct context_r *) context->custom_context;
	int n = context->n;
	double r = cr->r;

	return (n / 2.0) * log(2.0 * M_PI)
			+ (1.0 / (2.0 * r)) * cec_matrix_trace(cov) + (n / 2.0) * log(r);
}

double h_diagonal(const struct energy_function_context * context,
		const struct cec_matrix * cov)
{

	int n = context->n;

	double diagonal_product = cec_matrix_diagonal_product(cov);
	if (diagonal_product <= 0)
	{
		*context->last_error = POSITIVE_DEFINITE_ERROR;
		return NAN;
	}

	return (n / 2.0) * log(2.0 * M_PI * M_E)
			+ (1.0 / 2.0) * log(diagonal_product);
}

double h_all(const struct energy_function_context * context,
		const struct cec_matrix * cov)
{

	struct cec_matrix * temp_matrix = context->temp_matrix;
	int n = context->n;

	double det = cec_matrix_det_positive_definite(cov, temp_matrix);

	if (isnan(det))
	{
		*context->last_error = POSITIVE_DEFINITE_ERROR;
		return NAN;
	}
	return (n / 2.0) * log(2.0 * M_PI * M_E) + (1.0 / 2.0) * log(det);
}
