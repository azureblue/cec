#ifndef CEC_STARTER_H
#define CEC_STARTER_H

#include "vec.h"
#include "models/model.h"

namespace cec {

    using std::vector;
    using std::unique_ptr;
    using std::shared_ptr;

    class clustering_results;

    class cec_starter {
    public:
        cec_starter() = default;

        clustering_results start(const mat &x, const vector<int> &initial_assignment,
                                   const vector<unique_ptr<model>> &models, int max_iter, int min_card);

        vector<mat> split_points(const mat &points, const vector<int> &assignment, int k);
    };

    class clustering_results {
        friend class cec_starter;

    public:
        mat centers;
        vector<int> assignment;
        int cluster_number;
        int iterations;
        double energy;
        vector<mat> covariances;

        clustering_results(const clustering_results &res) = default;

    private:
        clustering_results(mat centers, vector<int> assignment, int cluster_number,
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