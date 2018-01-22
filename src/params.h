#ifndef CEC_PARAMS_H
#define CEC_PARAMS_H

#include "vec.h"
#include "models/model.h"
#include "models/all.h"
#include "models/spherical.h"
#include "models/diagonal.h"
#include "models/fixed_radius.h"
#include "models/covariance.h"
#include "models/eigenvalues.h"
#include "init.h"
#include "common.h"

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

    init_method parse_init_method(const string &method);

    model_type parse_model_type(const string &name);

    class centers_param {
    public:
        const init_method init_m;
        const mat centers_mat;
        const vector<int> var_centers;

        centers_param(init_method init_m, mat centers_mat, vector<int> var_centers)
                : init_m(init_m),
                  centers_mat(std::move(centers_mat)),
                  var_centers(std::move(var_centers)) {}

        shared_ptr<centers_init_spec> get_centers_init();
    };

    class control_param {
    public:
        int starts;
        int max_iter;
        int min_card;
        int threads;

        control_param(int starts, int max_iter, int min_card, int threads)
                : starts(starts),
                  max_iter(max_iter),
                  min_card(min_card),
                  threads(threads) {}
    };

    class split_param {
    public:
        int max_k;
        int max_depth;
        int tries;
        int initial_starts;

        split_param(int max_k, int max_depth, int tries, int initial_starts)
                : max_k(max_k),
                  max_depth(max_depth),
                  tries(tries),
                  initial_starts(initial_starts) {}
    };

    class model_spec {
    public:
        const model_type type;

        explicit model_spec(const model_type type)
                : type(type) {}

        virtual ~model_spec() = default;

        virtual unique_ptr<model> create_model() const = 0;

        static vector<unique_ptr<model>> create_models(vector<shared_ptr<model_spec>> specs);
        static vector<unique_ptr<model>> create_models(const model_spec &spec, int n = 1);
    };


    class model_all_spec : public model_spec {
    public:
        const int n;

        explicit model_all_spec(int n)
                : model_spec(model_type::ALL),
                  n(n) {}

        unique_ptr<model> create_model() const override {
            return make_unique<all>(n);
        }
    };

    class model_spherical_spec : public model_spec {
    public:
        const int n;

        explicit model_spherical_spec(int n)
                : model_spec(model_type::ALL),
                  n(n) {}

        unique_ptr<model> create_model() const override {
            return make_unique<spherical>(n);
        }
    };

    class model_diagonal_spec : public model_spec {
    public:
        const int n;

        explicit model_diagonal_spec(int n)
                : model_spec(model_type::ALL),
                  n(n) {}

        unique_ptr<model> create_model() const override {
            return make_unique<diagonal>(n);
        }
    };

    class model_fixed_radius_spec : public model_spec {
    public:
        const int n;
        const double r;

        explicit model_fixed_radius_spec(int n, double r)
                : model_spec(model_type::FIXED_R),
                  n(n),
                  r(r) {}

        unique_ptr<model> create_model() const override {
            return make_unique<fixed_radius>(n, r);
        }
    };

    class model_covariance_spec : public model_spec {
    public:
        const int n;
        const mat g_cov;

        explicit model_covariance_spec(int n, mat g_cov)
                : model_spec(model_type::FIXED_R),
                  n(n),
                  g_cov(std::move(g_cov)) {}

        unique_ptr<model> create_model() const override {
            return make_unique<covariance>(n, g_cov);
        }
    };

    class model_eigenvalues_spec : public model_spec {
    public:
        const int n;
        const vector<double> values;

        explicit model_eigenvalues_spec(int n, vector<double> values)
                : model_spec(model_type::EIGENVALUES),
                  n(n),
                  values(std::move(values)) {}

        unique_ptr<model> create_model() const override {
            return make_unique<eigenvalues>(n, values);
        }
    };

    class models_param {
    public:
        const vector<shared_ptr<model_spec>> specs;

        explicit models_param(vector<shared_ptr<model_spec>> specs)
                : specs(std::move(specs)) {}
    };
    
    
}
#endif //CEC_PARAMS_H
