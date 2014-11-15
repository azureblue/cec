#include <stdlib.h>

#include "alloc.h"
#include "errors.h"
#include "kmeanspp.h"
#include "rand.h"

static int binary_search_d(double key, double * array, int len)
{
    int a = 0, b = len - 1;
    while (a != b)
    {
	int mid = (a + b) / 2;
	if (array[mid] < key)
	{
	    a = mid + 1;
	} else
	{
	    b = mid;
	}
    }
    return a;
}

int kmeanspp(struct cec_matrix * X, struct cec_matrix * C)
{
    int m = X->m;
    int k = C->m;
    int n = X->n;

    double * dists = m_alloc(sizeof (double) * m);
    double * sums = m_alloc(sizeof (double) * m);

    if (dists == NULL || sums == NULL)
    {
	m_free(dists);
	m_free(sums);
	return MALLOC_ERROR;
    }

    array_copy(cec_matrix_row(X, 0), cec_matrix_row(C, 0), n);
    dists[0] = 0.0;
    sums[0] = 0.0;

    for (int i = 1; i < m; i++)
    {
	double dist = dist2(cec_matrix_row(X, i), cec_matrix_row(C, 0), n);
	dists[i] = dist;
	sums[i] = sums[i - 1] + dist;
    }
    
    cec_rand_init();

    for (int i = 1; i < k; i++)
    {
	double r = cec_rand();
	double n_sum = r * sums[m - 1];
	int idx = binary_search_d(n_sum, sums, m);
	array_copy(cec_matrix_row(X, idx), cec_matrix_row(C, i), n);
	
	for (int j = 1; j < m; j++)
	{
	    double dist = dist2(cec_matrix_row(X, j), cec_matrix_row(C, i), n);
	    if (dist < dists[j])
	    {
		dists[j] = dist;
	    }
	    sums[j] = sums[j - 1] + dists[j];
	}
    }
    
    cec_rand_end();
    
    m_free(dists);
    m_free(sums);
    return 0;
}
