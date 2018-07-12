#include "params.h"
#include "exceptions.h"
#include "models/all.h"
#include "models/spherical.h"
#include "models/diagonal.h"
#include "models/fixed_radius.h"
#include "models/fixed_covariance.h"
#include "models/fixed_eigenvalues.h"
#include "models/fixed_mean.h"

namespace cec {
    init_method parse_init_method(const string &method) {
        if (method == "none")
            return init_method::NONE;
        if (method == "kmeanspp")
            return init_method::KMEANSPP;
        if (method == "random")
            return init_method::RANDOM;
        throw invalid_init_method(method);
    }

    model_type parse_model_type(const string &name) {
        if (name == "all")
            return model_type::ALL;
        if (name == "covariance")
            return model_type::COVARIANCE;
        if (name == "diagonal")
            return model_type::DIAGONAL;
        if (name == "eigenvalues")
            return model_type::EIGENVALUES;
        if ((name == "radius") || (name == "fixed_r") || (name == "fixedr"))
            return model_type::FIXED_R;
        if (name == "mean")
            return model_type::MEAN;
        if (name == "spherical")
            return model_type::SPHERICAL;
        throw invalid_model_name(name);
    }

    vector<unique_ptr<model>>
    model_spec::create_models(vector<shared_ptr<model_spec>> specs) {
        int size = specs.size();
        vector<unique_ptr<model>> models(size);
        for (int i = 0; i < size; ++i)
            models[i] = specs[i]->create_model();
        return models;
    }

    vector<unique_ptr<model>> model_spec::create_models(const model_spec &spec, int n) {
        vector<unique_ptr<model>> models(n);
        for (int i = 0; i < n; ++i)
            models[i] = spec.create_model();
        return models;
    }

    shared_ptr<centers_init_spec> centers_param::get_centers_init() {
        switch (init_m) {
            case init_method::NONE:
                return make_shared<fixed_init_spec>(centers_mat);
            case init_method::RANDOM:
                return make_shared<random_init_spec>();
            default:
                return make_shared<kmeanspp_init_spec>();
        }
    }

    unique_ptr<model> model_all_spec::create_model() const {
        return make_unique<all>(n);
    }

    unique_ptr<model> model_spherical_spec::create_model() const {
        return make_unique<spherical>(n);
    }

    unique_ptr<model> model_diagonal_spec::create_model() const {
        return make_unique<diagonal>(n);
    }

    unique_ptr<model> model_fixed_radius_spec::create_model() const {
        return make_unique<fixed_radius>(n, r);
    }

    unique_ptr<model> model_covariance_spec::create_model() const {
        return make_unique<fixed_covariance>(n, g_cov);
    }

    unique_ptr<model> model_eigenvalues_spec::create_model() const {
        return make_unique<fixed_eigenvalues>(n, values);
    }

    unique_ptr<model> model_mean_spec::create_model() const {
        vec mean_vec(n);
        for (int i = 0; i < n; i++)
            mean_vec[i] = mean[i];
        return make_unique<fixed_mean>(n, mean_vec);
    }
}
