#include "alloc.h"
#include "kmeanspp.h"
#include "rand.h"

struct kmeanspp_ctx {
    double * dists_m;
    double * sums_m;
};

typedef struct kmeanspp_ctx kmeanspp_ctx;

static int binary_search_d(double key, const double * array, int len)
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

static void init(centers_init_ctx ctx, const cec_mat *x, cec_mat *c)
{
    int m = x->m;
    int k = c->m;
    int n = x->n;
    
    double * dists_m = ((kmeanspp_ctx*) ctx)->dists_m;
    double * sums_m = ((kmeanspp_ctx*) ctx)->sums_m;

    int first_center = (int) (cec_rand() * m);
    array_copy(cec_matrix_const_row(x, first_center), cec_matrix_row(c, 0), n);
    dists_m[0] = 0.0;
    sums_m[0] = 0.0;
    
    double * dists = dists_m + 1;
    double * sums = sums_m + 1;
    
    for (int i = 0; i < m; i++)
    {
        double dist = dist_sq(cec_matrix_const_row(x, i), cec_matrix_const_row(c, 0), n);
        dists[i] = dist;
        sums[i] = sums[i - 1] + dist;
    }
    for (int i = 1; i < k; i++)
    {
        double r = cec_rand();
        double n_sum = r * sums[m - 1];
        int idx = binary_search_d(n_sum, sums, m);
        array_copy(cec_matrix_const_row(x, idx), cec_matrix_row(c, i), n);

        for (int j = 0; j < m; j++)
        {
            double dist = dist_sq(cec_matrix_const_row(x, j), cec_matrix_const_row(c, i), n);
            if (dist < dists[j])
                dists[j] = dist;
            sums_m[j + 1] = sums_m[j] + dists[j];
        }
    }
}

centers_init * create_kmeanspp_initializer(int m) {
    centers_init * ci = alloc(centers_init);
    kmeanspp_ctx * km_ctx = alloc(kmeanspp_ctx);
    km_ctx->dists_m = alloc_n(double, m + 1);
    km_ctx->sums_m = alloc_n(double, m + 1);
    ci->ctx = km_ctx;
    ci->init = init;
    return ci;
}
