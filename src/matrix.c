#include "alloc.h"
#include "matrix.h"

void cec_matrix_set(cec_mat * m, double val)
{
    array_fill(m->data, val, m->m * m->n);
}

void cec_matrix_mul(cec_mat * m, double val)
{
    array_mul(m->data, val, m->m * m->n);
}

void cec_matrix_add(cec_mat *restrict m1, const cec_mat *restrict m2)
{
    array_add(m1->data, m2->data, m1->m * m1->n);
}

void cec_matrix_sub(cec_mat *restrict m1, const cec_mat *restrict m2)
{
    array_sub(m1->data, m2->data, m1->m * m1->n);
}

void cec_matrix_copy_data(const cec_mat *restrict from, cec_mat *restrict to)
{
    array_copy(from->data, to->data, from->m * from->n);
}

cec_mat * cec_matrix_create(int m, int n)
{
    cec_mat * mat = alloc_fam(cec_mat, double, m * n);
    mat->m = m;
    mat->n = n;
    return mat;
}

cec_mat * cec_matrix_create_copy(const cec_mat * m)
{
    cec_mat * rm = cec_matrix_create(m->m, m->n);
    if (!rm)
        return NULL;
    
    cec_matrix_copy_data(m, rm);
    return rm;
}

void cec_vector_outer_product(const double *restrict vec, cec_mat * output_matrix, int n)
{
    for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
            cec_matrix_set_element(output_matrix, j, k, vec[j] * vec[k]);
}
