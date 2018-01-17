#ifndef CEC_STARTER_H
#define CEC_STARTER_H

#include <utility>

#include "vec.h"
#include "common.h"
#include "models/model.h"

namespace cec {
    class clustering_results;

    struct starter_params {
        starter_params(int max_iter, int min_card)
                : max_iter(max_iter),
                  min_card(min_card) {}

        int max_iter;
        int min_card;
    };

    class points_split {
    public:
        points_split(mat points, vector<int> mapping)
                : points_(std::move(points)),
                  mapping_(std::move(mapping)) {}

        static vector<points_split> split_points(const mat &points, const vector<int> &assignment, int k);

        const mat &points() const {
            return points_;
        }

        const vector<int> &mapping() const {
            return mapping_;
        }

    private:
        mat points_;
        vector<int> mapping_;
    };

    class cec_starter {
    public:
        cec_starter() = default;

        unique_ptr<clustering_results>
        start(const mat &x, const vector<int> &initial_assignment,
              const vector<unique_ptr<model>> &models, const starter_params &params);
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

#endif //CEC_STARTER_H