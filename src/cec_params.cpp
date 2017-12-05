#include "cec_params.h"
#include "exceptions.h"

cec::centers_init_method cec::parse_init_method(const std::string &method) {
    if (method == "none")
        return centers_init_method::NONE;
    if (method == "kmeanspp")
        return centers_init_method::KMEANSPP;
    if (method == "random")
        return centers_init_method::RANDOM;
    throw invalid_init_method(method);
}

cec::model_name cec::parse_model_name(const std::string &name) {
    if (name == "all")
        return cec::model_name::ALL;
    if (name == "covariance")
        return cec::model_name::COVARIANCE;
    if (name == "diagonal")
        return cec::model_name::DIAGONAL;
    if (name == "eigenvalues")
        return cec::model_name::EIGENVALUES;
    if ((name == "fixed_r") || (name == "fixedr"))
        return cec::model_name::FIXED_R;
    if (name == "spherical")
        return cec::model_name::SPHERICAL;
    throw invalid_model_name(name);
}
