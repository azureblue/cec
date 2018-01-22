#include "r_params.h"
#include "r_utils.h"

namespace cec {
    namespace r {
        r_ext_ptr<centers_param> get_centers_param(SEXP centers_param_r) {
            r_wrapper r_par(centers_param_r);
            init_method im = parse_init_method(r_par["init.method"].get<const char *>());
            auto centers = (im == init_method::NONE) ? r_par["mat"].get<r_ext_ptr<mat>>() : make_r_ext<mat>(0, 0);
            auto var_centers = r_par["var.centers"].get<r_ext_ptr<vector<int>>>();
            return make_r_ext<centers_param>(im, *centers, *var_centers);
        }

        r_ext_ptr<control_param> get_control_param(SEXP control_param_r) {
            r_wrapper r_par(control_param_r);
            return make_r_ext<control_param>(
                    r_par["starts"].get<int>(),
                    r_par["max.iters"].get<int>(),
                    r_par["min.card"].get<int>(),
                    r_par["threads"].get<int>()
            );
        }

        r_ext_ptr<models_param> get_models_param(SEXP models_param_r, int n) {
            r_wrapper r_models(models_param_r);
            int len = r_models.size();
            auto specs = make_r_ext<vector<shared_ptr<model_spec>>>();
            for (int i = 0; i < len; i++) {
                r_wrapper model_r = r_models[i];
                model_type type = parse_model_type(model_r["type"].get<const char *>());
                r_wrapper params_r = model_r["params"];
                switch (type) {
                    case model_type::ALL:
                        specs->push_back(make_shared<model_all_spec>(n));
                        break;
                    case model_type::SPHERICAL:
                        specs->push_back(make_shared<model_spherical_spec>(n));
                        break;
                    case model_type::DIAGONAL:
                        specs->push_back(make_shared<model_diagonal_spec>(n));
                        break;
                    case model_type::FIXED_R:
                        specs->push_back(make_shared<model_fixed_radius_spec>(n, params_r["r"].get<double>()));
                        break;
                    case model_type::COVARIANCE: {
                        auto cov = params_r["cov"].get<r_ext_ptr<mat>>();
                        specs->push_back(make_shared<model_covariance_spec>(n, *cov));
                        break;
                    }
                    case model_type::EIGENVALUES: {
                        auto evals = params_r["eigenvalues"].get<r_ext_ptr<vector<double>>>();
                        specs->push_back(make_shared<model_eigenvalues_spec>(n, *evals));
                    }
                }
            }
            return make_r_ext<models_param>(std::move(*specs));
        }

        r_ext_ptr<split_param> get_split_param(SEXP split_param_r) {
            r_wrapper r_par(split_param_r);
            return make_r_ext<split_param> (
                    r_par["limit"].get<int>(),
                    r_par["depth"].get<int>(),
                    r_par["tries"].get<int>(),
                    r_par["initial.starts"].get<int>()
            );
        }
    }
}

