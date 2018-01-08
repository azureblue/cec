#ifndef CEC_STARTER_H
#define CEC_STARTER_H

#include "vec.h"
#include "models/model.h"

namespace cec {

    using std::vector;
    using std::unique_ptr;
    using std::shared_ptr;

    class single_start_results;

    class cec_starter {
    public:
        cec_starter() = default;

        single_start_results start(const mat &x, const vector<int> &initial_assignment,
                                   const vector<unique_ptr<model>> &models, int max_iter, int min_card);

        vector<mat> split_points(const mat &points, const vector<int> &assignment, int k);
    };

    class single_start_results {
        friend class cec_starter;

    public:
        const mat centers;
        const vector<int> assignment;
        const int cluster_number;
        const int iterations;
        const double energy;
        const vector<mat> covariances;

        single_start_results(const single_start_results &res) = default;

    private:
        single_start_results(mat centers, vector<int> assignment, int cluster_number,
                             int iterations, double energy, vector<mat> covariances)
                : centers(std::move(centers)),
                  assignment(std::move(assignment)),
                  cluster_number(cluster_number),
                  iterations(iterations),
                  energy(energy),
                  covariances(std::move(covariances)) {}
    };


}

#endif //CEC_STARTER_H