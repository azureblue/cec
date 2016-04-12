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
    int last_error;
    void * specific_context;
};

/*
 * Custom contexts for specific implementations.
 */

struct context_gc
{
    struct cec_matrix * given_cov;
    struct cec_matrix * i_given_cov;
    struct cec_matrix * temp_matrix;
};

struct context_fe
{
    double * evals_given;
    double * evals_temp;
    double given_evals_product;
    struct cec_matrix * workspace;
    struct cec_matrix * temp_matrix;
};

struct context_all
{
    struct cec_matrix * temp_matrix;
};

struct context_r
{
    double r;
};

struct cross_entropy_context * create_cross_entropy_context_all(int n);
struct cross_entropy_context * create_cross_entropy_context_spherical(int n);
struct cross_entropy_context * create_cross_entropy_context_diagonal(int n);
struct cross_entropy_context * create_cross_entropy_context_fixedr(int n, double r);
struct cross_entropy_context * create_cross_entropy_context_covariance(int n, const struct cec_matrix * cov, const struct cec_matrix * cov_i);
struct cross_entropy_context * create_cross_entropy_context_eigenvalues(int n, const double * given_evals);

void destroy_cross_entropy_context_all(struct cross_entropy_context * context);
void destroy_cross_entropy_context_spherical(struct cross_entropy_context * context);
void destroy_cross_entropy_context_diagonal(struct cross_entropy_context * context);
void destroy_cross_entropy_context_fixedr(struct cross_entropy_context * context);
void destroy_cross_entropy_context_covariance(struct cross_entropy_context * context);
void destroy_cross_entropy_context_eigenvalues(struct cross_entropy_context * context);

#endif	/* CROSS_ENTROPY_CONTEXT_H */
