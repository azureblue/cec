#ifndef ENERGY_FUNCTION_CONTEXT_H
#define	ENERGY_FUNCTION_CONTEXT_H

enum density_family
{
    GIVEN_COVARIANCE = 0, FIXED_R = 1, SPHERICAL = 2, DIAGONAL = 3, FIXEDEIGENVALUES = 4, ALL = 5
};

struct energy_function_context
{
    int n;
    struct cec_matrix * temp_matrix;
    int * last_error;
    void * custom_context;
    enum density_family family;
};

struct context_gc
{
    struct cec_matrix * given_cov;
    struct cec_matrix * i_given_cov;
};

struct context_fe
{
    double * given_evals;
    double given_evals_product;
    double * evals;
    struct cec_matrix * workspace;
};

struct context_r
{
    double r;
};

void destroy_energy_function_context(struct energy_function_context * context);

struct energy_function_context * create_energy_function_context(
	enum density_family family, int n);

#endif	/* ENERGY_FUNCTION_CONTEXT_H */

