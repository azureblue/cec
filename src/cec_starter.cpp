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

    std::vector<int> assignment = in.initial_assignment;

    std::vector<mat> split = split_points(in.x, in.initial_assignment, k);
    int removed_k = 0;

    double energy_sum = 0;
    int iterations = 0;

    std::vector<std::unique_ptr<cluster>> clusters(k);

    for (int i = 0; i < k; i++) {
        mat &cluster_split = split[i];
        if (cluster_split.m >= min_card) {
            mean me(cluster_split);
            mat cov = covariance_mle::estimate(cluster_split, me);
            clusters[i].reset(new cluster(*in.models[i], me, cov, m));
        }
    }

    for (int i = 0; i < k; i++) {
        if (!clusters[i])
            continue;

        double energy = clusters[i]->energy();

        if (std::isnan(energy))
            throw invalid_covariance(*clusters[i], i);

        energy_sum += energy;
    }

    int removed = std::count(clusters.begin(), clusters.end(), std::unique_ptr<cluster>());
    if (removed == k)
        throw all_clusters_removed();

    bool handle_removed_flag = removed_k != 0;

    for (int iter = (handle_removed_flag ? -1 : 0); iter < max_iter; iter++) {
        bool transfer_flag = false;
        bool removed_last_iteration_flag = false;

        for (int i = 0; i < m; i++) {

            const int cl_num = assignment[i];
            std::unique_ptr<cluster> &cl_src = clusters[cl_num];

            if (handle_removed_flag && cl_src)
                continue;

            double rem_energy_gain = cl_src ? cl_src->rem_point(X[i]) : 0;
            double best_gain = cl_src ? 0 : std::numeric_limits<double>::infinity();

            if (std::isnan(rem_energy_gain))
                throw invalid_covariance(*cl_src, cl_num);

            int dest_cl_num = -1;

            for (int j = 0; j < k; j++) {
                if (j == cl_num || !clusters[j])
                    continue;

                cluster &cl_dst = *clusters[j];
                double add_energy_gain = cl_dst.add_point(X[i]);

                if (std::isnan(add_energy_gain))
                    throw invalid_covariance(cl_dst, j);

                double gain = add_energy_gain + rem_energy_gain;

                if (gain < best_gain) {
                    dest_cl_num = j;
                    best_gain = gain;
                }

            }

            if (dest_cl_num != -1) {
                assignment[i] = dest_cl_num;
                clusters[dest_cl_num]->apply_change();
                if (cl_src) {
                    cl_src->apply_change();
                    if (cl_src->card() < min_card) {
                        removed_last_iteration_flag = true;
                        energy_sum -= cl_src->energy();
                        clusters[cl_num].reset(nullptr);
                    }
                }
                energy_sum += best_gain;
                transfer_flag = true;
            }
        }

        iterations = iter + 1;
        if (!transfer_flag)
            break;

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
        if (!clusters[i]) {
            centers[i].fill(std::numeric_limits<double>::quiet_NaN());
            covs[i].fill(std::numeric_limits<double>::quiet_NaN());
            continue;
        }
        centers[i] = clusters[i]->mean();
        covs[i] = clusters[i]->covariance();
        energy_check += clusters[i]->energy();
    }
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

