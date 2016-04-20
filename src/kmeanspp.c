#include "alloc.h"
#include "errors.h"
#include "kmeanspp.h"
#include "rand.h"

static int binary_search_d(double key, double * array, int len)
{
    int a = 0, b = len - 1;
    while (a != b)
    {
        if (array[a] == 0 && array[b] == 0)
            return (int) (cec_rand() * (b - a + 1)) + a;
        int mid = (a + b) / 2;
        if (array[mid] < key)
            a = mid + 1;        
        else
            b = mid;
    }
    return a;
}

int kmeanspp(struct cec_matrix * X, struct cec_matrix * C)
{
    int m = X->m;
    int k = C->m;
    int n = X->n;
    
    double * dists_m = m_alloc(sizeof (double) * (m + 1));
    double * sums_m = m_alloc(sizeof (double) * (m + 1));
    
    if (!dists_m || !sums_m)
    {
        m_free(dists_m);
        m_free(sums_m);
        return MALLOC_ERROR;
    }
    
    cec_rand_init();
    
    int first_center = (int) (cec_rand() * m);
    
    array_copy(cec_matrix_const_row(X, first_center), cec_matrix_row(C, 0), n);
    dists_m[0] = 0.0;
    sums_m[0] = 0.0;
    
    double * dists = dists_m + 1;
    double * sums = sums_m + 1;
    
    for (int i = 0; i < m; i++)
    {
        double dist = dist2(cec_matrix_const_row(X, i), cec_matrix_const_row(C, 0), n);
        dists[i] = dist;
        sums[i] = sums[i - 1] + dist;
    }

    for (int i = 1; i < k; i++)
    {
        double r = cec_rand();
        double n_sum = r * sums[m - 1];
        int idx = binary_search_d(n_sum, sums, m);
        array_copy(cec_matrix_const_row(X, idx), cec_matrix_row(C, i), n);

        for (int j = 0; j < m; j++)
        {
            double dist = dist2(cec_matrix_const_row(X, j), cec_matrix_const_row(C, i), n);
            if (dist < dists[j])
                dists[j] = dist;
            
            sums_m[j + 1] = sums_m[j] + dists[j];
        }
    }
    
    cec_rand_end();
    
    m_free(dists_m);
    m_free(sums_m);
    return 0;
}
