#include "split_starter.h"
#include "exceptions.h"

namespace cec {
    unique_ptr<clustering_results>
    split_starter::try_split_cluster(const mat &x_mat) {
        try {
            const unique_ptr<clustering_results> &single_res
                    = cec.start(x_mat, vector<int>(x_mat.m, 0),
                                model_spec::create_models(m_spec));

            if (!single_res)
                return unique_ptr<clustering_results>();

            unique_ptr<clustering_results> split_res = splitter.start(clustering_input(x_mat, try_split_models));

            if (split_res && split_res->cluster_number == 2 &&
                split_res->energy < single_res->energy)
                return split_res;

            return unique_ptr<clustering_results>();
        } catch (clustering_exception &ce) {
            return unique_ptr<clustering_results>();
        }
    }

    unique_ptr<clustering_results>
    split_starter::start(const unique_ptr<clustering_results> &cl_res,
                         const clustering_input &input_params) {

        if (cl_res->cluster_number >= max_k)
            return make_unique<clustering_results>(*cl_res);

        const mat &x = input_params.x;
        int n = x.n;
        int k = cl_res->centers.m;
        vector<int> cluster(x.m);
        vector<bool> moved(k, true);

        auto current_res = make_unique<clustering_results>(*cl_res);

        for (int split_level = 0; split_level < max_depth; split_level++) {
            cluster = current_res->assignment;
            k = current_res->centers.m;
            const vector<points_split> &split = points_split::split_points(x, cluster, k);
            mat split_res_c_mat = mat(2, n);
            int k_s = 0;
            int next_k = std::min(max_k, k * 2);
            vector<bool> is_split(next_k);
            bool split_flag = false;

            for (int i = 0; i < k; i++) {
                if (k_s >= max_k)
                    break;
                int split_m = split[i].points().m;
                if (split_m == 0)
                    continue;
                bool split_success = false;

                const mat &split_x_mat = split[i].points();
                const vector<int> mapping = split[i].mapping();
                vector<int> split_assignment(split_x_mat.m);
                if (moved[i]) {
                    auto const &split_res = try_split_cluster(split_x_mat);
                    if (split_res) {
                        split_success = true;
                        split_res_c_mat = split_res->centers;
                        split_assignment = split_res->assignment;
                    }
                }
                split_flag = split_flag || split_success;
                if (k_s == max_k - 1)
                    split_success = false;

                if (split_success) {
                    is_split[k_s] = true;
                    is_split[k_s + 1] = true;

                    for (int p = 0; p < split_m; p++) {
                        if (split_assignment[p] == 0)
                            cluster[mapping[p]] = k_s;
                        else
                            cluster[mapping[p]] = k_s + 1;
                    }
                    k_s += 2;

                } else {
                    is_split[k_s] = false;
                    for (int p = 0; p < split_m; p++)
                        cluster[mapping[p]] = k_s;
                    k_s++;
                }
            }

            bool need_another_split = false;
            if (split_flag) {
                try {
                    current_res = cec.start(x, cluster, model_spec::create_models(m_spec, k_s));
                } catch (clustering_exception &ce) {
                    return current_res;
                }

                for (int i = 0; i < x.m; i++) {
                    if (cluster[i] != current_res->assignment[i]) {
                        moved[cluster[i]] = true;
                        moved[current_res->assignment[i]] = true;
                    }
                }
                for (int i = 0; i < k_s; i++) {
                    if (is_split[i])
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
