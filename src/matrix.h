#ifndef MATRIX_DOUBLE_H
#define	MATRIX_DOUBLE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "errors.h"
#include "array.h"

struct cec_matrix
{
    int m;
    int n;
    double * data;
};

struct cec_matrix * cec_matrix_create(int m, int n);

void cec_matrix_destroy(struct cec_matrix *restrict m);

void cec_matrix_set(struct cec_matrix *restrict m, double val);

static inline void cec_matrix_set_element(struct cec_matrix *restrict m, int a,
	int b, double val)
{
    m->data[m->n * a + b] = val;
}

static inline double cec_matrix_element(const struct cec_matrix *restrict m,
	int a, int b)
{
    return m->data[m->n * a + b];
}

static inline double * cec_matrix_row(const struct cec_matrix *restrict matrix,
	int m)
{
    return matrix->data + (m * matrix->n);
}

void cec_matrix_mul(struct cec_matrix *restrict m, double val);

void cec_matrix_add(struct cec_matrix *restrict m1,
	const struct cec_matrix *restrict m2);

void cec_matrix_sub(struct cec_matrix *restrict m1,
	const struct cec_matrix *restrict m2);

void cec_matrix_sum_multiplied(const struct cec_matrix *restrict m1,
	const struct cec_matrix *restrict m2, struct cec_matrix * dest,
	double a1, double a2);

void cec_matrix_copy_data(const struct cec_matrix *restrict from,
	struct cec_matrix *restrict to);

void cec_vector_outer_product(const double *restrict vec,
	struct cec_matrix *restrict output_matrix, int n);


#endif	/* MATRIX_DOUBLE_H */

