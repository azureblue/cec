#include "cec_r.h"
#include "r_utils.h"
#include "params.h"
#include "r_params.h"
#include "starter.h"
#include "r_result.h"
#include "cec_starter.h"
#include "variable_starter.h"
#include "parallel_starter.h"

#include<R_ext/Random.h>

using namespace cec;
using namespace cec::r;
using std::exception;

static void seed_from_r() {
    GetRNGstate();
    double r = unif_rand();
    PutRNGstate();
    unsigned long seed;
    memcpy(&seed, &r, sizeof(seed));
    random::set_seed(seed);
}

#define ERR_STR_LEN 1000
extern "C"
SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r) {

    seed_from_r();
    bool exception_thrown = false;
    r_ext_ptr<std::array<char, ERR_STR_LEN + 1>> error_str_array;
    try {
        error_str_array.init();
    } catch (exception &ex) {
        exception_thrown = true;
    }

    if (exception_thrown) {
        error("cec library error");
    }

    r_ext_ptr<clustering_results> start_results;
    try {
        mat x_mat = get<mat>(x);
        int n = x_mat.n;

        centers_param centers_par = get_centers_param(centers_param_r);
        control_param control_par = get_control_param(control_param_r);
        models_param models_par = get_models_param(models_param_r, n);

        const shared_ptr<centers_init_spec> &centers_init_ptr = centers_par.get_centers_init();
        const centers_init_spec &init_spec = *centers_init_ptr;

        //clustering_input input_params(x_mat, model_spec::create_models(models_par.specs));

        cec_parameters start_params(control_par.max_iter, control_par.min_card);
        cec_starter::parameters initial_cl_params(start_params, init_spec);

        split_starter::parameters split_params(start_params, *models_par.specs[0], init_spec, 10,
                                               20, 10);
        start_and_split_task split_task(initial_cl_params, split_params);

        multiple_starts_task m_start_task(initial_cl_params);

        parallel_starter ps(control_par.threads, control_par.starts);
        unique_ptr<clustering_results> results = ps.start(split_task, x_mat, models_par.specs);

        start_results.reset(results.release());

    } catch (exception &ex) {
        exception_thrown = true;
        std::strncat(error_str_array->data(), ex.what(), ERR_STR_LEN);
    }

    if (exception_thrown) {
        error(error_str_array->data());
    }

    try {
        SEXP r_res;
        PROTECT(r_res = create_R_result(*start_results));
        UNPROTECT(1);
        return r_res;
    } catch (std::exception &ex) {
        std::strncat(error_str_array->data(), ex.what(), ERR_STR_LEN);
    }
    error(error_str_array->data());
}

extern "C"
SEXP cec_init_centers_r(SEXP x_r, SEXP k_r, SEXP method_r) {
    seed_from_r();
    r_ext_ptr<std::array<char, ERR_STR_LEN>> error_str_array;
    error_str_array.init();
    r_ext_ptr<mat> res;
    try {
        try {
            const string &method_str = get<string>(method_r);
            const mat &x = get<mat>(x_r);
            int k = get<int>(k_r);
            init_method im = parse_init_method(method_str);
            switch (im) {
                case init_method::KMEANSPP:
                    res.init(kmeanspp_init().init(x, k));
                    break;
                case init_method::RANDOM:
                    res.init(random_init().init(x, k));
                    break;
                default:
                    throw invalid_init_method(method_str);
            }
        } catch (std::exception &ex) {
            throw;
        }
        SEXP r_res = put(*res);
        return r_res;
    } catch (exception &ex) {
        std::strncat(error_str_array->data(), ex.what(), ERR_STR_LEN);
    }

    error(error_str_array->data());
}

R_CallMethodDef methods[] = {
        {"cec_r",              (DL_FUNC) &cec_r,              4},
        {"cec_init_centers_r", (DL_FUNC) &cec_init_centers_r, 3},
        {NULL, NULL,                                          0}
};

extern "C"
void R_init_CEC(DllInfo *dllInfo) {
    R_registerRoutines(dllInfo, NULL, methods, NULL, NULL);
    R_useDynamicSymbols(dllInfo, TRUE);
}
