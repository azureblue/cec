#ifndef MATRIX_H
#define	MATRIX_H

#include "errors.h"
#include "array.h"

struct cec_matrix
{
    int m;
    int n;
    double data[];
};

typedef struct cec_matrix cec_mat;

cec_mat * cec_matrix_create(int m, int n);

cec_mat * cec_matrix_create_copy(const cec_mat * mat);

void cec_matrix_destroy(cec_mat * m);

void cec_matrix_set(cec_mat * m, double val);

static inline void cec_matrix_set_element(cec_mat * m, int a, int b, double val)
{
    m->data[m->n * a + b] = val;
}

static inline double cec_matrix_element(const cec_mat * m, int a, int b)
{
    return m->data[m->n * a + b];
}

static inline const double * cec_matrix_const_row(const cec_mat * matrix, int m)
{
    return matrix->data + (m * matrix->n);
}

static inline double * cec_matrix_row(cec_mat * matrix, int m)
{
    return matrix->data + (m * matrix->n);
}

void cec_matrix_mul(cec_mat * m, double val);

void cec_matrix_add(cec_mat *restrict m1, const cec_mat *restrict m2);

void cec_matrix_sub(cec_mat *restrict m1, const cec_mat *restrict m2);

void cec_matrix_copy_data(const cec_mat *restrict from, cec_mat *restrict to);

void cec_vector_outer_product(const double *restrict vec, cec_mat *restrict output_matrix, int n);

#endif	/* MATRIX_H */
