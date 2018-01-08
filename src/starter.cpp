#include "starter.h"
#include "cluster.h"
#include "exceptions.h"
#include "models/constants.h"

cec::single_start_results cec::cec_starter::start(const mat &x, const vector<int> &initial_assignment,
                                                  const vector<unique_ptr<model>> &models, int max_iter, int min_card) {

    int m = x.m;
    int k = models.size();
    int n = x.n;
    vector<int> assignment = initial_assignment;
    double energy_sum = 0;
    int iters = 0;

    vector<unique_ptr<cluster>> clusters(k);
    vector<mat> split = split_points(x, assignment, k);

    for (int i = 0; i < k; i++) {
        mat &cluster_split = split[i];
        if (cluster_split.m >= min_card) {
            mean me(cluster_split);
            mat cov = covariance_mle::estimate(cluster_split, me);
            clusters[i].reset(new cluster(*models[i], me, cov, m));
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

    int removed_after_assignmenet = std::count(clusters.begin(), clusters.end(), unique_ptr<cluster>());
    if (removed_after_assignmenet == k)
        throw all_clusters_removed();

    bool handle_removed_flag = removed_after_assignmenet != 0;

    for (int iter = (handle_removed_flag ? -1 : 0); iter < max_iter; iter++) {
        bool transfer_flag = false;
        bool removed_last_iteration_flag = false;

        for (int i = 0; i < m; i++) {

            const int cl_num = assignment[i];
            unique_ptr<cluster> &cl_src = clusters[cl_num];

            if (handle_removed_flag && cl_src)
                continue;

            double rem_energy_gain = cl_src ? cl_src->rem_point(x[i]) : 0;
            double best_gain = cl_src ? 0 : constants::INF;

            if (std::isnan(rem_energy_gain))
                throw invalid_covariance(*cl_src, cl_num);

            int dest_cl_num = -1;

            for (int j = 0; j < k; j++) {
                if (j == cl_num || !clusters[j])
                    continue;

                cluster &cl_dst = *clusters[j];
                double add_energy_gain = cl_dst.add_point(x[i]);

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

        iters = iter + 1;
        if (!transfer_flag)
            break;

        if (removed_last_iteration_flag) {
            handle_removed_flag = true;
            iter--;
        } else
            handle_removed_flag = false;
    }

    mat centers(k, n);
    vector<mat> covs(k, mat(n, n));

    for (int i = 0; i < k; i++) {
        if (!clusters[i]) {
            centers[i].fill(constants::QNAN);
            covs[i].fill(constants::QNAN);
            continue;
        }
        centers[i] = clusters[i]->mean();
        covs[i] = clusters[i]->covariance();
    }
    int final_k = k - std::count(clusters.begin(), clusters.end(), unique_ptr<cluster>());
    return single_start_results(centers, assignment, final_k, iters, energy_sum, covs);
}

std::vector<cec::mat>
cec::cec_starter::split_points(const mat &points, const vector<int> &assignment, int k) {
    vector<int> sizes(k, 0);
    vector<int> indices(k, 0);
    for (auto &&cl : assignment)
        sizes[cl]++;
    vector<mat> split;
    for (auto &&size : sizes)
        split.emplace_back(size, points.n);
    int m = points.m;
    for (int i = 0; i < m; i++) {
        int cl = assignment[i];
        split[cl][indices[cl]++] = points[i];
    }
    return split;
}

