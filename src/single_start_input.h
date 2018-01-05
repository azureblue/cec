#ifndef CEC_SINGLE_START_INPUT_H
#define CEC_SINGLE_START_INPUT_H

#include <utility>

#include "vec.h"
#include "models/all.h"

namespace cec {
    class cec_starter;
    class single_start_results {
    friend class cec_starter;
    public:
        const mat centers;
        const std::vector<int> assignment;
        const int cluster_number;
        const int iterations;
        const double energy;
        const std::vector<mat> covariances;
        single_start_results(const single_start_results &res) = default;

    private:
        single_start_results(mat centers, std::vector<int> assignment, int cluster_number,
                             int iterations, double energy, std::vector<mat> covariances)
                : centers(std::move(centers)),
                  assignment(std::move(assignment)),
                  cluster_number(cluster_number),
                  iterations(iterations),
                  energy(energy),
                  covariances(std::move(covariances)) {}
    };

    class single_start_input {
    public:
        const mat &x, &c;
        const std::vector<int> &initial_assignment;
        const std::vector<std::unique_ptr<model>> &models;
        const int max_iter;
        const int min_card;

        single_start_input(const mat &x, const mat &c, const std::vector<int> &initial_assignment,
                           const std::vector<std::unique_ptr<model>> &models, const int max_iter, const int min_card) :
                x(x), c(c),
                initial_assignment(initial_assignment),
                models(models),
                max_iter(max_iter),
                min_card(min_card) {}
    };
}
#endif //CEC_SINGLE_START_INPUT_H
