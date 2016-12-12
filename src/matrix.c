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

void cec_matrix_add(struct cec_matrix *restrict m1, const struct cec_matrix *restrict m2)
{
    array_add(m1->data, m2->data, m1->m * m1->n);
}

void cec_matrix_sub(struct cec_matrix *restrict m1, const struct cec_matrix *restrict m2)
{
    array_sub(m1->data, m2->data, m1->m * m1->n);
}

void cec_matrix_sum_multiplied(struct cec_matrix *restrict dest, 
        const struct cec_matrix *restrict m1, double a1, const struct cec_matrix *restrict m2, double a2)
{
    array_sum_multiplied(dest->data, m1->data, a1, m2->data, a2, m1->m * m1->n);
}

void cec_matrix_copy_data(const struct cec_matrix *restrict from, struct cec_matrix *restrict to)
{
    array_copy(from->data, to->data, from->m * from->n);
}

struct cec_matrix * cec_matrix_create(int m, int n)
{
    struct cec_matrix * mat = alloc_fam(struct cec_matrix, double, m * n);
    
    if (!mat)
        return NULL;    
    
    mat->m = m;
    mat->n = n;
    return mat;
}

struct cec_matrix * cec_matrix_create_copy(const struct cec_matrix * m)
{
    struct cec_matrix * rm = cec_matrix_create(m->m, m->n);
    if (!rm)
        return NULL;
    
    cec_matrix_copy_data(m, rm);
    return rm;
}

void cec_matrix_destroy(struct cec_matrix * m)
{
    m_free(m);
}

void cec_vector_outer_product(const double *restrict vec, struct cec_matrix * output_matrix, int n)
{
    for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
            cec_matrix_set_element(output_matrix, j, k, vec[j] * vec[k]);
}
