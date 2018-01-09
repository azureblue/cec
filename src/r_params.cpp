#include "r_params.h"
#include "r_utils.h"

using namespace cec::r;
using std::shared_ptr;

cec::centers_param cec::get_centers_param(SEXP centers_param_r) {
    r_wrapper r_par(centers_param_r);
    init_method im = parse_init_method(r_par["init.method"].get<string>());
    auto centers = r_ext_ptr<mat>::make((im == init_method::NONE)
                                        ? r_par["mat"].get<mat>()
                                        : mat(0, 0));
    const vector<int> &var_centers = r_par["var.centers"].get<vector<int>>();
    return centers_param(im, *centers, var_centers);
}

cec::control_param cec::get_control_param(SEXP control_param_r) {
    r_wrapper r_par(control_param_r);
    return {
            r_par["starts"].get<int>(),
            r_par["max.iters"].get<int>(),
            r_par["min.card"].get<int>(),
            0
    };
}

cec::models_param cec::get_models_param(SEXP models_param_r, int n) {
    r_wrapper r_models(models_param_r);
    int len = r_models.size();
    vector<std::shared_ptr<model_spec>> specs;
    for (int i = 0; i < len; i++) {
        r_wrapper model_r = r_models[i];
        const string &type_name = model_r["type"].get<string>();
        model_type type = parse_model_type(type_name);
        r_wrapper params_r = model_r["params"];
        switch (type) {
            case model_type::ALL:
                specs.push_back(shared_ptr<model_spec>(new model_all_spec(n)));
                break;
            case model_type::SPHERICAL:
                specs.push_back(shared_ptr<model_spec>(new model_spherical_spec(n)));
                break;
            case model_type::DIAGONAL:
                specs.push_back(shared_ptr<model_spec>(new model_diagonal_spec(n)));
                break;
            case model_type::FIXED_R:
                specs.push_back(shared_ptr<model_spec>(
                        new model_fixed_radius_spec(n, params_r["r"].get<double>())));
                break;
            case model_type::COVARIANCE:
                specs.push_back(shared_ptr<model_spec>(
                        new model_covariance_spec(n, params_r["cov"].get<mat>())));
                break;
            case model_type::EIGENVALUES:
                specs.push_back(shared_ptr<model_spec>(
                        new model_eigenvalues_spec(n, params_r["eigenvalues"].get<vector<double>>())
                ));
        }
    }
    return models_param(std::move(specs));
}

