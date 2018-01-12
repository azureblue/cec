#include "params.h"
#include "exceptions.h"

cec::init_method cec::parse_init_method(const string &method) {
    if (method == "none")
        return init_method::NONE;
    if (method == "kmeanspp")
        return init_method::KMEANSPP;
    if (method == "random")
        return init_method::RANDOM;
    throw invalid_init_method(method);
}

cec::model_type cec::parse_model_type(const string &name) {
    if (name == "all")
        return cec::model_type::ALL;
    if (name == "covariance")
        return cec::model_type::COVARIANCE;
    if (name == "diagonal")
        return cec::model_type::DIAGONAL;
    if (name == "eigenvalues")
        return cec::model_type::EIGENVALUES;
    if ((name == "fixed_r") || (name == "fixedr"))
        return cec::model_type::FIXED_R;
    if (name == "spherical")
        return cec::model_type::SPHERICAL;
    throw invalid_model_name(name);
}

std::vector<std::unique_ptr<cec::model>>
cec::model_spec::create_models(std::vector<std::shared_ptr<cec::model_spec>> specs) {
    int size = specs.size();
    vector<unique_ptr<model>> models(size);
    for (int i = 0; i < size; ++i)
        models[i] = specs[i]->create_model();
    return models;
}

std::shared_ptr<cec::centers_init_spec> cec::centers_param::get_centers_init() {
    switch (init_m) {
        case init_method::NONE:
            return shared_ptr<centers_init_spec>(new fixed_init_spec(centers_mat));
        case init_method::RANDOM:
            return shared_ptr<centers_init_spec>(new random_init_spec());
        default:
            return shared_ptr<centers_init_spec>(new kmeanspp_init_spec());
    }
}
