#include "cec_r.h"
#include "r_utils.h"
#include "params.h"
#include "r_params.h"
#include "init.h"
#include "starter.h"
#include "r_result.h"
#include "multi_starter.h"
#include<R_ext/Random.h>

using namespace cec;
using namespace cec::r;
using std::vector;
using std::string;
using std::unique_ptr;
using std::shared_ptr;

static void seed_r() {
    GetRNGstate();
    double r = unif_rand();
    PutRNGstate();
    unsigned long seed;
    memcpy(&seed, &r, sizeof(seed));
    random::context::set_seed(seed);
}

extern "C"
SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r) {
    seed_r();
    r_ext_ptr<clustering_results> start_results;
    try {
        mat x_mat = get<mat>(x);
        int n = x_mat.n;

        centers_param centers_par = get_centers_param(centers_param_r);
        control_param control_par = get_control_param(control_param_r);
        models_param models_par = get_models_param(models_param_r, n);

        int k = centers_par.var_centers[0];
        vector<unique_ptr<model>> models(k);

        std::transform(models_par.specs.begin(), models_par.specs.end(), models.begin(),
                       [&](const shared_ptr<cec::model_spec> &spec) {
                           return spec->create_model();
                       });
        multi_starter ms;
        init_spec is(shared_ptr<centers_init_spec>(new kmeanspp_init_spec(x_mat.m)),
                     shared_ptr<assignment_init_spec>(new closest_init_spec()));

        unique_ptr<clustering_results> results = ms.start(x_mat, models_par.specs, is, control_par.max_iterations,
                                                          control_par.min_card, control_par.starts,
                                                          control_par.threads);
        start_results.reset(results.release());

    } catch (std::exception &ex) {
        error(ex.what());
    }

    try {
        SEXP r_res;
        PROTECT(r_res = create_R_result(*start_results));
        UNPROTECT(1);
        return r_res;
    } catch (std::exception &ex) {
        error(ex.what());
    }
}

extern "C"
SEXP cec_init_centers_r(SEXP x_r, SEXP k_r, SEXP method_r) {
    seed_r();
    r_ext_ptr<string> error_str;
    try {
        r_ext_ptr<mat> res;
        try {
            const string &method_str = get<string>(method_r);
            const mat &x = get<mat>(x_r);
            int k = get<int>(k_r);
            init_method im = parse_init_method(method_str);
            switch (im) {
                case init_method::KMEANSPP:
                    res.reset(new mat(kmeanspp_init(x.m).init(x, k)));
                    break;
                case init_method::RANDOM:
                    res.reset(new mat(random_init().init(x, k)));
                    break;
                case init_method::NONE:
                    break;
            }
        } catch (std::exception &ex) {
            throw;
        }
        SEXP r_res = put(*res);
        return r_res;
    } catch (std::exception &ex) {
        error_str.reset(new string(ex.what()));
    }
    error(error_str->c_str());
}
