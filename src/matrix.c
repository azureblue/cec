#include "alloc.h"
#include "matrix.h"

void cec_matrix_set(struct cec_matrix * m, double val)
{
	array_fill(m->data, val, m->m * m->n);
}

void cec_matrix_mul(struct cec_matrix * m, double val)
{
	array_mul(m->data, val, m->m * m->n);
}

void cec_matrix_add(struct cec_matrix * m1, const struct cec_matrix * m2)
{
	array_add(m1->data, m2->data, m1->m * m1->n);
}

void cec_matrix_sub(struct cec_matrix * m1, const struct cec_matrix * m2)
{
	array_sub(m1->data, m2->data, m1->m * m1->n);
}

void cec_matrix_sum_multiplied(const struct cec_matrix * m1,
		const struct cec_matrix * m2, struct cec_matrix * dest, double a1,
		double a2)
{
	array_sum_multiplied(dest->data, m1->data, a1, m2->data, a2, m1->m * m1->n);
}

void cec_matrix_copy_data(const struct cec_matrix * from,
		struct cec_matrix * to)
{
	array_copy(from->data, to->data, from->m * from->n);
}

struct cec_matrix * cec_matrix_create(int m, int n)
{
	struct cec_matrix * mat = (struct cec_matrix *) m_alloc(
			sizeof(struct cec_matrix) + sizeof(double) * m * n);
	if (mat == NULL)
		return NULL;
	mat->data = (double *) (mat + 1);
	mat->m = m;
	mat->n = n;
	return mat;
}

void cec_matrix_destroy(struct cec_matrix * m)
{
	m_free(m);
}

void cec_vector_outer_product(const double * vec,
		struct cec_matrix * output_matrix, int n)
{
	for (int j = 0; j < n; j++)
		for (int k = 0; k < n; k++)
			cec_matrix_set_element(output_matrix, j, k, vec[j] * vec[k]);
}

double cec_matrix_trace(const struct cec_matrix * m)
{
	double res = 0;
	for (int i = 0; i < m->n; i++)
	{
		res += cec_matrix_element(m, i, i);
	}
	return res;
}

double cec_matrix_trace_assert_positive(const struct cec_matrix * m)
{
	double res = 0;
	for (int i = 0; i < m->n; i++)
	{
		double el = cec_matrix_element(m, i, i);
		if (el <= 0)
		{
			return NAN;
		}
		res += el;
	}
	return res;
}

void cec_matrix_multiply_sq(const struct cec_matrix * m1,
		const struct cec_matrix * m2, struct cec_matrix * dest)
{
	/*
	 Simple cec_matrix multiplication.
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

int cec_cholesky(const struct cec_matrix * m, struct cec_matrix * target)
{
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

double cec_matrix_diagonal_product(const struct cec_matrix * m)
{
	double res = 1.0;
	int n = m->n;
	for (int i = 0; i < n; i++)
		res *= cec_matrix_element(m, i, i);
	return res;
}

double cec_matrix_det_positive_definite(const struct cec_matrix * m,
		struct cec_matrix * temp)
{
	/*
	 * Special case for 2x2 matrix.
	 */

	if (m->n == 2)
	{
		double det = m->data[0] * m->data[3] - m->data[1] * m->data[2];
		if (det > 0 && m->data[0] + m->data[3])
			return det;
		else
			return NAN;
	} else if (cec_cholesky(m, temp) == POSITIVE_DEFINITE_ERROR)
	{
		return NAN;
	}
	double prod = cec_matrix_diagonal_product(temp);
	return prod * prod;
}

