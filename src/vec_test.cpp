#include "vec.h"
#include "cluster.h"
#include "cec_starter.h"

using namespace cec;

int main() {


    cec_starter starter(2);
    int m = 70000;
    int n = 2;
    int k = 3;
    mat pts(m, n);
    mat cen(k, n);
    for (int i = 0; i < m; i++) {
        pts[i][0] = rand() / (double) RAND_MAX;
        pts[i][1] = rand() / (double) RAND_MAX;
    }

    for (int i = 0; i < k; i++) {
        cen[i] = {rand() / (double) RAND_MAX, rand() / (double) RAND_MAX};
    }

    std::vector<int> asgn(m);
    for (int i = 0; i < m; i++) {
        asgn[i] = rand() % k;
    }

    std::vector<std::unique_ptr<model>> models;
    models.push_back(std::unique_ptr<model>(new all(n)));
    models.push_back(std::unique_ptr<model>(new all(n)));
    models.push_back(std::unique_ptr<model>(new all(n)));
    single_start_input in(pts, cen, asgn, models, 10, 10);
    auto results = starter.start(in);
    std::cout << results.energy << std::endl;
    std::cout << results.iterations << std::endl;
    for (auto &&cl : results.assignment) {
        std::cout << cl << " ";
    }
}