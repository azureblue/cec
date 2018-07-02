#include "init.h"

std::vector<int> cec::closest_assignment::init(const cec::mat &x, const cec::mat &c) {
    int m = x.m;
    int k = c.m;
    std::vector<int> asgn(m);
    for (int i = 0; i < m; i++) {
        double b_dist = std::numeric_limits<double>::max();
        const row &point = x[i];
        int b_row = -1;
        for (int j = 0; j < k; j++) {
            double dist = row::dist(point, c[j]);
            if (dist < b_dist) {
                b_dist = dist;
                b_row = j;
            }
        }
        asgn[i] = b_row;
    }
    return asgn;
}

cec::mat cec::random_init::init(const mat &x, int k) {
    mat c_mat(k, x.n);
    std::uniform_int_distribution<int> unif_int(0, x.m - 1);
    for (int i = 0; i < k; i++)
        c_mat[i] = x[unif_int(mt)];
    return c_mat;
}

cec::mat cec::kmeanspp_init::init(const mat &x, int k) {
    int m = x.m;
    int n = x.n;
    mat c(k, n);
    std::uniform_int_distribution<int> unif_int(0, x.m - 1);
    c[0] = x[unif_int(mt)];
    dists[0] = 0.0;
    sums[0] = 0.0;

    sums[0] = dists[0] = row::dist(x[0], c[0]);

    for (int i = 1; i < m; i++) {
        double dist = row::dist(x[i], c[0]);
        dists[i] = dist;
        sums[i] = sums[i - 1] + dist;
    }

    std::uniform_real_distribution<double> unif_real;
    for (int i = 1; i < k; i++) {
        double upper = sums[m - 1];
        double n_sum = upper == 0.0
                       ? 0.0
                       : unif_real(mt, std::uniform_real_distribution<double>::param_type(0.0, upper));
        auto range = std::equal_range(sums.begin(), sums.end(), n_sum);
        int idx_from = range.first - sums.begin();
        int idx_to = range.second - sums.begin();
        idx_to = std::min(idx_to, m - 1);
        int idx = unif_int(mt, std::uniform_int_distribution<int>::param_type(idx_from, idx_to));
        c[i] = x[idx];
        sums[0] = dists[0] = std::min(dists[0], row::dist(x[0], c[i]));
        for (int j = 1; j < m; j++) {
            dists[j] = std::min(dists[j], row::dist(x[j], c[i]));
            sums[j] = sums[j - 1] + dists[j];
        }
    }
    return c;
}
