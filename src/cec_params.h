#ifndef CEC_PARAMS_H
#define CEC_PARAMS_H

#include "vec.h"
#include "models/model.h"
#include "models/all.h"
#include "models/spherical.h"

namespace cec {
    enum class init_method {
        NONE,
        KMEANSPP,
        RANDOM
    };

    enum class model_type {
        ALL,
        COVARIANCE,
        DIAGONAL,
        EIGENVALUES,
        FIXED_R,
        SPHERICAL
    };

    init_method parse_init_method(const std::string &method);

    model_type parse_model_type(const std::string &name);

    class centers_param {
    public:
        const init_method init_m;
        const mat centers_mat;
        const std::vector<int> var_centers;

        centers_param(init_method init_m, mat centers_mat, std::vector<int> var_centers)
                : init_m(init_m),
                  centers_mat(std::move(centers_mat)),
                  var_centers(std::move(var_centers)) {}
    };

    class control_param {
    public:
        const int starts;
        const int max_iterations;
        const int min_card;
        const int threads;

        control_param(int starts, int max_iter, int min_card, int threads)
                : starts(starts),
                  max_iterations(max_iter),
                  min_card(min_card),
                  threads(threads) {}
    };

    class model_spec {
    public:
        const model_type type;

        virtual std::unique_ptr<model> create_model() const = 0;

        explicit model_spec(const model_type type)
                : type(type) {}
    };

    class model_all_spec : public model_spec {
    public:
        const int n;

        explicit model_all_spec(int n)
                : model_spec(model_type::ALL),
                  n(n) {}

        std::unique_ptr<model> create_model() const override {
            return std::unique_ptr<model>(new all(n));
        }
    };

    class model_spherical_spec : public model_spec {
    public:
        const int n;

        explicit model_spherical_spec(int n)
                : model_spec(model_type::ALL),
                  n(n) {}

        std::unique_ptr<model> create_model() const override {
            return std::unique_ptr<model>(new spherical(n));
        }
    };

    class models_param {
    public:
        const std::vector<std::shared_ptr<model_spec>> specs;

        explicit models_param(std::vector<std::shared_ptr<model_spec>> specs)
                : specs(std::move(specs)) {}
    };
}
#endif //CEC_PARAMS_H
