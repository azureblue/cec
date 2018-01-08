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
