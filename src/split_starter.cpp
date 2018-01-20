#include "split_starter.h"

namespace cec {
    unique_ptr<clustering_results>
    split_starter::try_split_cluster(const mat &x_mat, shared_ptr<model_spec> m_spec) {

        const unique_ptr<clustering_results> &single_res
                = cec.start(x_mat, vector<int>(x_mat.m, 0), model_spec::create_models(m_spec));

        if (!single_res)
            return unique_ptr<clustering_results>();

        double single_cluster_energy = single_res->energy;

        unique_ptr<clustering_results> split_res
                = splitter.start(x_mat, vector<shared_ptr<model_spec>>({m_spec, m_spec}));

        if (split_res && split_res->cluster_number == 2 &&
            split_res->energy < single_cluster_energy)
            return split_res;

        return unique_ptr<clustering_results>();
    }

    unique_ptr<clustering_results>
    split_starter::start(const unique_ptr<clustering_results> &cl_res, const mat &x_mat,
                         const shared_ptr<model_spec> &model_sp) {

        int n = x_mat.n;
        int k = cl_res->centers.m;
        int split_depth = params.max_depth;
        vector<int> cluster(x_mat.m);
        vector<bool> moved(k, true);

        auto current_res = make_unique<clustering_results>(*cl_res);

        for (int split_level = 0; split_level < split_depth; split_level++) {
            cluster = current_res->assignment;
            k = current_res->centers.m;
            const vector<points_split> &split = points_split::split_points(x_mat, cluster, k);
            mat split_centers = mat(k * 2, n);
            mat split_res_c_mat = mat(2, n);
            int k_s = 0;
            vector<bool> splitted(k * 2);
            bool split_flag = false;

            for (int i = 0; i < k; i++) {
                int split_m = split[i].points().m;
                if (split_m == 0)
                    continue;
                split_centers[k_s] = current_res->centers[i];
                bool split_success = false;

                const mat &split_x_mat = split[i].points();
                const vector<int> mapping = split[i].mapping();
                vector<int> split_assignment(split_x_mat.m);
                if (moved[i]) {
                    const unique_ptr<clustering_results> &split_res = try_split_cluster(
                            split_x_mat,
                            model_sp);
                    if (split_res) {
                        split_success = true;
                        split_res_c_mat = split_res->centers;
                        split_assignment = split_res->assignment;
                    }
                }
                split_flag = split_flag || split_success;
                if (split_success) {
                    splitted[k_s] = true;
                    splitted[k_s + 1] = true;
                    split_centers[k_s] = split_res_c_mat[0];
                    split_centers[k_s + 1] = split_res_c_mat[1];

                    for (int p = 0; p < split_m; p++) {
                        if (split_assignment[p] == 0)
                            cluster[mapping[p]] = k_s;
                        else
                            cluster[mapping[p]] = k_s + 1;
                    }
                    k_s += 2;

                } else {
                    splitted[k_s] = false;
                    split_centers[k_s] = current_res->centers[i];
                    for (int p = 0; p < split_m; p++)
                        cluster[mapping[p]] = k_s;
                    k_s++;
                }
            }

            mat level_c_mat(k_s, n);
            std::copy_n(split_centers.begin(), k_s, level_c_mat.begin());
            bool need_another_split = false;
            if (split_flag) {
                current_res = cec.start(x_mat, cluster, model_spec::create_models(model_sp));

                for (int i = 0; i < x_mat.m; i++) {
                    if (cluster[i] != current_res->assignment[i]) {
                        moved[cluster[i]] = true;
                        moved[current_res->assignment[i]] = true;
                    }
                }
                for (int i = 0; i < k_s; i++) {
                    if (splitted[i])
                        moved[i] = true;
                    need_another_split = (need_another_split || moved[i]);
                }
            }
            if (!need_another_split)
                break;
        }
        return current_res;
    }
}