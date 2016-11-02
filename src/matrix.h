#ifndef MATRIX_H
#define	MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "errors.h"
#include "array.h"

/*
 * Simple matrix structure used in CEC and basic set of matrix functions.
 */

struct cec_matrix
{
    int m;
    int n;
    double data[];
};

struct cec_matrix * cec_matrix_create(int m, int n);

struct cec_matrix * cec_matrix_create_copy(const struct cec_matrix * mat);

void cec_matrix_destroy(struct cec_matrix * m);

void cec_matrix_set(struct cec_matrix * m, double val);

static inline void cec_matrix_set_element(struct cec_matrix * m, int a, int b, double val)
{
    m->data[m->n * a + b] = val;
}

static inline double cec_matrix_element(const struct cec_matrix * m, int a, int b)
{
    return m->data[m->n * a + b];
}

static inline const double * cec_matrix_const_row(const struct cec_matrix * matrix, int m)
{
    return matrix->data + (m * matrix->n);
}

static inline double * cec_matrix_row(struct cec_matrix * matrix, int m)
{
    return matrix->data + (m * matrix->n);
}

void cec_matrix_mul(struct cec_matrix * m, double val);

void cec_matrix_add(struct cec_matrix *restrict m1, const struct cec_matrix *restrict m2);

void cec_matrix_sub(struct cec_matrix *restrict m1, const struct cec_matrix *restrict m2);

void cec_matrix_sum_multiplied(struct cec_matrix * dest, 
        const struct cec_matrix *restrict m1, double a1, 
        const struct cec_matrix *restrict m2, double a2);

void cec_matrix_copy_data(const struct cec_matrix *restrict from, struct cec_matrix *restrict to);

void cec_vector_outer_product(const double *restrict vec, struct cec_matrix *restrict output_matrix, int n);

#endif	/* MATRIX_H */
