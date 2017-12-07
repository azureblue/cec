#include "cec_params_r.h"
#include "r_utils.h"

using cec::r::r_wrapper;

cec::centers_param cec::get_centers_param(SEXP centers_param_r) {
    r_wrapper r_par = centers_param_r;
    init_method im = parse_init_method(r_par.get<std::string>("init.method"));
    mat centers = (im == init_method::NONE)
                  ? r_par.get<mat>("mat")
                  : mat(0, 0);
    centers_param cent_param(im, centers, r_par.get<std::vector<int>>("var.centers"));
}

cec::control_param  cec::get_control_param(SEXP control_param_r) {
    r_wrapper r_par = control_param_r;
    return {
            r_par.get<int>("starts"),
            r_par.get<int>("max.iters"),
            r_par.get<int>("min.card"),
            r_par.get<int>("threads")
    };
}

cec::models_param cec::get_models_param(SEXP models_param_r, int n) {
    r_wrapper r_models = models_param_r;
    int len = r_models.size();
    std::vector<std::shared_ptr<model_spec>> specs;
    for (int i = 0; i < len; i++) {
        r_wrapper model_r = r_models.get<SEXP>(i);
        const std::string &type_name = model_r.get<std::string>("type");
        model_type type = parse_model_type(type_name);
        r_wrapper params_r = r_models.get<SEXP>("params");
        switch (type) {
            case model_type::ALL:
                specs.push_back(std::shared_ptr<model_spec>(
                        new model_all_spec(n)
                ));
                break;
            case model_type::COVARIANCE:
            case model_type::DIAGONAL:
            case model_type::EIGENVALUES:
            case model_type::FIXED_R:
            case model_type::SPHERICAL:
                throw not_implemented(type_name);
        }
    }
    return models_param(std::move(specs));
}

