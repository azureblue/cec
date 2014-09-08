#include "covariance.h"

void covariance_add_point(const struct cec_matrix * covariance,
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

void covariance_remove_point(const struct cec_matrix * covariance,
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

