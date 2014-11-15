#ifndef CROSS_ENTROPY_CONTEXT_H
#define	CROSS_ENTROPY_CONTEXT_H

enum density_family
{
    GIVEN_COVARIANCE = 0, FIXED_R = 1, SPHERICAL = 2, DIAGONAL = 3, FIXEDEIGENVALUES = 4, ALL = 5
};

/*
 * Static context of all cross entropy functions.
 */
struct cross_entropy_context
{
    int n;
    struct cec_matrix * temp_matrix;
    int * last_error;
    void * custom_context;
    enum density_family family;
};

/*
 * Custom contexts for specific implementations.
 */

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

void destroy_cross_entropy_context(struct cross_entropy_context * context);

struct cross_entropy_context * create_cross_entropy_context(
	enum density_family family, int n);

#endif	/* CROSS_ENTROPY_CONTEXT_H */

