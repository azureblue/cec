#ifndef CEC_PARAMS_H
#define CEC_PARAMS_H

#include "vec.h"

namespace cec {
    enum class centers_init_method {
        NONE,
        KMEANSPP,
        RANDOM,
    };

    enum class model_name {
        ALL, COVARIANCE, DIAGONAL, EIGENVALUES, FIXED_R, SPHERICAL
    };

    centers_init_method parse_init_method(const std::string &method);

    model_name parse_model_name(const std::string &name);

    struct cec_centers_param {
        centers_init_method init_m;
        mat centers_mat;
        std::vector<int> var_centers;

        cec_centers_param(centers_init_method init_m, mat centers_mat, std::vector<int> var_centers):
                init_m(init_m),
                centers_mat(std::move(centers_mat)),
                var_centers(std::move(var_centers)) {}

    };

    struct cec_control_param {
        int starts;
        int max_iterations;
        int min_card;
        int threads;

        cec_control_param(int starts, int max_iterations, int min_card, int threads):
                starts(starts),
                max_iterations(max_iterations),
                min_card(min_card),
                threads(threads) {}
    };

    struct cec_model_r_params {
        double r;
    };

    struct cec_model_eigenvalues_params {
        const std::vector<double> *given_eigenvalues;
    };

    struct cec_model_covariances_params {
        const mat cov;
        const mat cov_inv;
    };

    struct cec_model_specification {
        enum model_name type;
        int n;
        union {
            struct cec_model_r_params r_params;
            struct cec_model_eigenvalues_params eigenvalues_params;
            struct cec_model_covariances_params covariances_params;
        };
    };

    struct cec_models_param {
        int len;
        struct cec_model_specification *model_specs;
    };
}
#endif //CEC_PARAMS_H
