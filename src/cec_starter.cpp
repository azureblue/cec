#include <map>
#include <iomanip>
#include <cmath>
#include "cec_starter.h"
#include "cluster.h"
#include "exceptions.h"

cec::single_start_results cec::cec_starter::start(const cec::single_start_input &in) {

    std::cout << std::setprecision(20);
    const mat &X = in.x;
    int m = X.m;
    int k = in.c.m;
    int n = X.n;
    int max_iter = in.max_iter;
    int min_card = in.min_card;

    const std::unique_ptr<model> *models = in.models.data();
    std::vector<int> assignment = in.initial_assignment;

    std::vector<mat> split = split_points(in.x, in.initial_assignment, k);
    std::vector<bool> removed(k, false);
    int removed_k = 0;

    double energy_sum = 0;
    int iterations = 0;

    std::vector<std::unique_ptr<cluster>> clusters(k);

    for (int i = 0; i < k; i++) {
        const mat &cluster_split = split[i];
        if (cluster_split.m < min_card) {
            removed_k++;
            removed[i] = true;
            continue;
        }
        mean me(cluster_split);
        mat cov = covariance_mle::estimate(cluster_split, me);
        clusters[i].reset(new cluster(*models[i], me, cov, m));
    }

    for (int i = 0; i < k; i++) {
        if (removed[i])
            continue;

        double energy = clusters[i]->energy();

        if (std::isnan(energy))
            throw invalid_covariance(*clusters[i], i);

        energy_sum += energy;
    }

    std::cout << energy_sum << std::endl;

    for (auto &cl : clusters) {

        std::cout << cl->mean() << std::endl;
        std::cout << cl->covariance() << std::endl;
    }

    if (removed_k == k)
        throw all_clusters_removed();

    bool handle_removed_flag = removed_k != 0;

    for (int iter = (handle_removed_flag ? -1 : 0); iter < max_iter; iter++) {
        int transfer_flag = 0;
        int removed_last_iteration_flag = 0;

        for (int i = 0; i < m; i++) {

            const int cl_num = assignment[i];
            const bool cl_removed = removed[cl_num];
            cluster &cl = *clusters[cl_num];

            if (handle_removed_flag && !cl_removed)
                continue;

            double rem_energy_gain = !cl_removed ? cl.rem_point(X[i]) : 0;
            double best_gain = !cl_removed ? 0 : std::numeric_limits<double>::infinity();
            if (std::isnan(rem_energy_gain))
                throw invalid_covariance(cl, cl_num);

            int dest_cl_num = -1;

            for (int j = 0; j < k; j++) {
                if ((j == cl_num) || removed[j])
                    continue;

                cluster &cl_j = *clusters[j];
                double add_energy_gain = cl_j.add_point(X[i]);

                if (std::isnan(add_energy_gain))
                    throw invalid_covariance(cl_j, j);

                double gain = add_energy_gain;
                if (!cl_removed)
                    gain += rem_energy_gain;

                if (gain < best_gain) {
                    dest_cl_num = j;
                    best_gain = gain;
                }

            }

            if (dest_cl_num != -1) {
                assignment[i] = dest_cl_num;
                clusters[dest_cl_num]->apply_change();
                if (!cl_removed) {
                    cl.apply_change();
                    if (cl.card() < min_card) {
                        removed_last_iteration_flag = 1;

                        energy_sum -= cl.energy();
                        removed[cl_num] = true;
                    }
                }
                energy_sum += best_gain;
                transfer_flag = 1;
            }
        }

        iterations = iter + 1;
        std::cout << energy_sum << std::endl;
        if (!transfer_flag)
            break;

        /*
         * If cluster was removed in this iteration, we need to perform another iteration
         * considering points that are not assigned. It will prevent energy dips.
         */

        if (removed_last_iteration_flag) {
            handle_removed_flag = true;
            iter--;
        } else
            handle_removed_flag = false;
    }

    mat centers = in.c;
    std::vector<mat> covs(centers.m, mat(n, n));

    double energy_check = 0;
    for (int i = 0; i < k; i++) {
        if (removed[i])
            continue;
        centers[i] = clusters[i]->mean();
        covs[i] = clusters[i]->covariance();
        energy_check += clusters[i]->energy();
    }
    std::cout << energy_check << " - check" << std::endl;
    return single_start_results(centers, assignment, k, iterations, energy_sum, covs);

}

std::vector<cec::mat>
cec::cec_starter::split_points(const mat &points, const std::vector<int> &assignment, int k) {
    std::vector<int> sizes(k, 0);
    std::vector<int> indices(k, 0);
    for (auto &&cl : assignment)
        sizes[cl]++;
    std::vector<mat> split;
    for (auto &&size : sizes)
        split.emplace_back(size, points.n);
    int m = points.m;
    for (int i = 0; i < m; i++) {
        int cl = assignment[i];
        split[cl][indices[cl]++] = points[i];
    }
    return split;
}

