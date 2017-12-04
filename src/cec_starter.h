#ifndef CEC_STARTER_H
#define CEC_STARTER_H

#include "vec.h"
#include "single_start_input.h"

namespace cec {

    class cec_starter {
    public:
        cec_starter(int n): t_mean_matrix(n, n), t_matrix_nn(n, n), n_covariance_matrix(n, n) {}
        single_start_results start(const single_start_input &in);
        std::vector<std::vector<vec>> split_points(const mat &points, const std::vector<int> &assignment);


    private:
        mat t_mean_matrix;
        mat t_matrix_nn;
        mat n_covariance_matrix;
        std::vector<mat> t_covariance_matrices;

        class cluster_map {
        public:
            const int k;

            explicit cluster_map(int k):
                    k(k),
                    removed_(k, false),
                    mapping_(k) {
                for (int i = 0; i < k; i++)
                    mapping_[i] = i;
            }

        private:
            int current_k_ = k;
            int removed_clusters_ = 0;
            std::vector<bool> removed_;
            std::vector<int> mapping_;
        };
    };


}

#endif //CEC_STARTER_H