#include "vec.h"
#include "cluster.h"
#include "cec_starter.h"

using namespace cec;

int main() {

    srand(1234);
    int m = 4500;
    int k = 3;
    int n = 3;
    mat x(m, n);
    FILE * xf = fopen("/media/docs_ssd/Projects/cec/inst/cec_tests/mouse3d.data", "r");
    if (!xf)
        exit(-1);
    for (int i = 0; i < 4500; ++i) {
        double xx, yy, zz;
        if (!fscanf(xf, "%lf %lf %lf", &xx, &yy, &zz))
            exit(-1);
        x[i][0] = xx;
        x[i][1] = yy;
        x[i][2] = zz;
    }
    fclose(xf);

    mat c(k, n);
    FILE * cf = fopen("/media/docs_ssd/Projects/cec/inst/cec_tests/centers3d.data", "r");
    if (!cf)
        exit(-1);
    for (int i = 0; i < k; i++) {
        double xx, yy, zz;
        if (!fscanf(xf, "%lf %lf %lf", &xx, &yy, &zz))
            exit(-1);
        c[i][0] = xx;
        c[i][1] = yy;
        c[i][2] = zz;
    }
    fclose(cf);

    std::vector<int> asgn(m);
    for (int i = 0; i < m; i++) {
        asgn[i] = i % k;
    }

    std::vector<std::unique_ptr<model>> models;
    models.push_back(std::unique_ptr<model>(new all(n)));
    models.push_back(std::unique_ptr<model>(new all(n)));
    models.push_back(std::unique_ptr<model>(new all(n)));
    single_start_input in(x, c, asgn, models, 50, 10);
    cec_starter starter(n);
    auto results = starter.start(in);
    std::cout << results.energy << std::endl;
    std::cout << results.iterations << std::endl;
    for (auto &&cl : results.assignment) {
        std::cout << cl << " ";
    }

    return 0;


//    cec_starter starter(2);
//    int m = 70000;
//    int n = 2;
//    int k = 3;
//    mat pts(m, n);
//    mat cen(k, n);
//    for (int i = 0; i < m; i++) {
//        pts[i][0] = rand() / (double) RAND_MAX;
//        pts[i][1] = rand() / (double) RAND_MAX;
//    }
//
//    for (int i = 0; i < k; i++) {
//        cen[i] = {rand() / (double) RAND_MAX, rand() / (double) RAND_MAX};
//    }
//
//    std::vector<int> asgn(m);
//    for (int i = 0; i < m; i++) {
//        asgn[i] = rand() % k;
//    }
//
//    std::vector<std::unique_ptr<model>> models;
//    models.push_back(std::unique_ptr<model>(new all(n)));
//    models.push_back(std::unique_ptr<model>(new all(n)));
//    models.push_back(std::unique_ptr<model>(new all(n)));
//    single_start_input in(pts, cen, asgn, models, 10, 10);
//    auto results = starter.start(in);
//    std::cout << results.energy << std::endl;
//    std::cout << results.iterations << std::endl;
//    for (auto &&cl : results.assignment) {
//        std::cout << cl << " ";
//    }
}