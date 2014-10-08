#include <stdio.h>
#include <float.h>
#include <R_ext/Lapack.h>

#include "cov_utils.h"

double cec_cov_trace(const struct cec_matrix * m)
{
    double res = 0;
    for (int i = 0; i < m->n; i++)
    {
	res += cec_matrix_element(m, i, i);
    }
    return res;
}

double cec_cov_diagonal_product(const struct cec_matrix * m)
{
    double res = 1.0;
    int n = m->n;
    for (int i = 0; i < n; i++)
	res *= cec_matrix_element(m, i, i);
    return res;
}

void cec_cov_multiply(const struct cec_matrix * m1,
	const struct cec_matrix * m2, struct cec_matrix * dest)
{
    /*
     Simple matrix multiplication.
     */
    int n = m1->n;
    cec_matrix_set(dest, 0.0);
    for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	    for (int k = 0; k < n; k++)
	    {
		cec_matrix_set_element(dest, i, j,
			cec_matrix_element(dest, i, j)
			+ cec_matrix_element(m1, k, i)
			* cec_matrix_element(m2, j, k));
	    }
}

int cec_cov_eigenvalues(const struct cec_matrix * sym_matrix, struct cec_matrix * temp_matrix, struct cec_matrix * workspace, double * values)
{
    int n = sym_matrix->n;
    int info;
    array_copy(sym_matrix->data, temp_matrix->data, n * n);
    int workspace_len = workspace->m * workspace->n;
    F77_NAME(dsyev)("N", "U", &n, temp_matrix->data, &n, values, workspace->data, &workspace_len, &info);
    if (info != 0)
	return UNKNOWN_ERROR;
    return 0;
}

int cec_cov_cholesky(const struct cec_matrix * sym_matrix, struct cec_matrix * temp_matrix){
    int n = sym_matrix->n;
    int info;
    array_copy(sym_matrix->data, temp_matrix->data, n * n);   
    F77_NAME(dpotrf)("U", &n, temp_matrix->data, &n, &info);  
    if (info != 0) 
	return POSITIVE_DEFINITE_ERROR;
    else
	return 0;
}
/* -- using lapack cholesky instead --
 
int cec_cov_cholesky(const struct cec_matrix * m, struct cec_matrix * target)
{
    return cec_cov_cholesky(m, target);
    int n = m->n;
    for (int i = 0; i < n; i++)
    {
	double sjk = 0.0;

	for (int j = 0; j < i; j++)
	{
	    double s = 0.0;

	    for (int k = 0; k < j; k++)
	    {
		s += cec_matrix_element(target, i, k)
			* cec_matrix_element(target, j, k);
	    }
	    double lij = (cec_matrix_element(m, i, j) - s)
		    / cec_matrix_element(target, j, j);
	    if (isinf(lij))
	    {
		return POSITIVE_DEFINITE_ERROR;
	    }

	    sjk += lij * lij;
	    cec_matrix_set_element(target, i, j, lij);
	}

	for (int j = i + 1; j < n; j++)
	{
	    cec_matrix_set_element(target, i, j, 0.0);
	}

	double lii = sqrt(cec_matrix_element(m, i, i) - sjk);

	if (isnan(lii))
	{
	    return POSITIVE_DEFINITE_ERROR;
	}

	cec_matrix_set_element(target, i, i, lii);
    }
    return NO_ERROR;
}
*/

double cec_cov_cholesky_det(const struct cec_matrix * m,
	struct cec_matrix * temp)
{
    /*
     * Special case for 2x2 matrix.
     */

    if (m->n == 2)
    {
	double det = m->data[0] * m->data[3] - m->data[1] * m->data[2];
	return det;
	
    } else if (cec_cov_cholesky(m, temp) == POSITIVE_DEFINITE_ERROR)
    {
	return NAN;
    }
    double prod = cec_cov_diagonal_product(temp);
    return prod * prod;
}

void cec_cov_add_point(const struct cec_matrix * covariance,
	struct cec_matrix * new_covarioance, const double * mean,
	double const * point, int card, struct cec_matrix * t_matrix)
{
    int n = covariance->n;
    double vec[n];
    array_copy(mean, vec, n);
    array_sub(vec, point, n);
    cec_vector_outer_product(vec, t_matrix, n);
    double d_card_1 = (double) card + 1.0;
    cec_matrix_sum_multiplied(covariance, t_matrix, new_covarioance,
	    card / d_card_1, card / (d_card_1 * d_card_1));
}

void cec_cov_remove_point(const struct cec_matrix * covariance,
	struct cec_matrix * new_covarioance, const double * mean,
	double const * point, int card, struct cec_matrix * t_matrix)
{
    int n = covariance->n;
    double vec[n];
    array_copy(mean, vec, n);
    array_sub(vec, point, n);
    cec_vector_outer_product(vec, t_matrix, n);
    double d_card_1 = (double) card - 1.0;
    cec_matrix_sum_multiplied(covariance, t_matrix, new_covarioance,
	    card / d_card_1, -card / (d_card_1 * d_card_1));
}
