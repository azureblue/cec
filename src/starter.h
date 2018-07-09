#ifndef CROSS_ENTROPY_CLUSTERING_H
#define CROSS_ENTROPY_CLUSTERING_H

#include <utility>

#include "vec.h"
#include "common.h"
#include "models/model.h"

namespace cec {
    class clustering_results;

    struct cec_parameters {
        cec_parameters(int max_iter, int min_card)
                : max_iter(max_iter),
                  min_card(min_card) {}

        int max_iter;
        int min_card;
    };

    class points_split {
    public:
        static vector<points_split> split_points(const mat &points, const vector<int> &assignment, int k);

        const mat &points() const {
            return pts;
        }

        const vector<int> &mapping() const {
            return map;
        }

        points_split(mat points, vector<int> mapping)
                : pts(std::move(points)),
                  map(std::move(mapping)) {}

    private:
        mat pts;
        vector<int> map;
    };

    class cross_entropy_clustering {
    public:
        explicit cross_entropy_clustering(const cec_parameters &params)
                : params(params) {}

        unique_ptr<clustering_results>
        start(const mat &x, const vector<int> &initial_assignment,
              const vector<unique_ptr<model>> &models);

    private:
        cec_parameters params;
    };

    class clustering_results {
        friend class cross_entropy_clustering;

    public:
        mat centers;
        vector<int> assignment;
        int cluster_number;
        int iterations;
        double energy;
        vector<mat> covariances;

        clustering_results(const clustering_results &res) = default;

        clustering_results(clustering_results &&res) = default;

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

#endif //CROSS_ENTROPY_CLUSTERING_H
