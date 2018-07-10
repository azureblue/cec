#include "params.h"
#include "starter.h"
#include "cec_starter.h"
#include "variable_starter.h"

#include "r_result.h"
#include "r_utils.h"
#include "r_params.h"
#include "cec_r.h"

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

extern "C"
SEXP cec_r(SEXP x_r, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r) {
    seed_from_r();
    const char *ex_what = nullptr;
    r_ext_ptr<clustering_results> start_results;
    try {
        auto x = get<r_ext_ptr<mat>>(x_r);
        int n = x->n;

        r_ext_ptr<centers_param> centers_par = get_centers_param(centers_param_r);
        r_ext_ptr<control_param> control_par = get_control_param(control_param_r);
        r_ext_ptr<models_param> models_par = get_models_param(models_param_r, n);

        const shared_ptr<centers_init_spec> &centers_init_ptr = centers_par->get_centers_init();
        const centers_init_spec &init_spec = *centers_init_ptr;

        cec_starter::parameters starter_params({control_par->max_iter, control_par->min_card}, init_spec);

        parallel_starter ps(control_par->threads, control_par->starts);

        auto cl_function = [&](const mat &x, const vector<shared_ptr<model_spec>> &specs) {
            multiple_starts_task task(starter_params, x, specs);
            return ps.start(task);
        };

        variable_starter var_start(std::move(cl_function), centers_par->var_centers);

        unique_ptr<clustering_results> results = var_start.start(*x, models_par->specs);

        if (!results)
            throw clustering_exception("all starts failed");

        start_results.reset(results.release());

    } catch (exception &ex) {
        ex_what = ex.what();
    }

    if (ex_what)
        error(ex_what);

    try {
        SEXP r_res;
        PROTECT(r_res = create_R_result(*start_results));
        UNPROTECT(1);
        return r_res;
    } catch (std::exception &ex) {
        ex_what = ex.what();
    }
    error(ex_what);
}

extern "C"
SEXP cec_split_r(SEXP x_r, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r,
                 SEXP split_param_r) {
    seed_from_r();
    const char *ex_what = nullptr;
    r_ext_ptr<clustering_results> start_results;
    try {
        auto x = get<r_ext_ptr<mat>>(x_r);
        int n = x->n;

        r_ext_ptr<centers_param> centers_par = get_centers_param(centers_param_r);
        r_ext_ptr<control_param> control_par = get_control_param(control_param_r);
        r_ext_ptr<models_param> models_par = get_models_param(models_param_r, n);
        r_ext_ptr<split_param> split_par = get_split_param(split_param_r);

        const shared_ptr<centers_init_spec> &centers_init_ptr = centers_par->get_centers_init();
        const centers_init_spec &init_spec = *centers_init_ptr;

        cec_parameters start_params(control_par->max_iter, control_par->min_card);

        cec_starter::parameters initial_cl_params(start_params, init_spec,
                                                  split_par->initial_starts);

        split_starter::parameters split_params(start_params, *models_par->specs[0], init_spec,
                                               split_par->tries, split_par->max_k,
                                               split_par->max_depth);

        parallel_starter ps(control_par->threads, control_par->starts);

        auto cl_function = [&](const mat &x, const vector<shared_ptr<model_spec>> &specs) {
            start_and_split_task task(initial_cl_params, split_params, x, specs);
            return ps.start(task);
        };

        variable_starter var_start(cl_function, centers_par->var_centers);

        unique_ptr<clustering_results> results = var_start.start(*x, models_par->specs);

        start_results.reset(results.release());

    } catch (exception &ex) {
        ex_what = ex.what();
    }

    if (ex_what)
        error(ex_what);

    try {
        SEXP r_res;
        PROTECT(r_res = create_R_result(*start_results));
        UNPROTECT(1);
        return r_res;
    } catch (std::exception &ex) {
        ex_what = ex.what();
    }
    error(ex_what);
}

extern "C"
SEXP cec_init_centers_r(SEXP x_r, SEXP k_r, SEXP method_r) {
    seed_from_r();
    const char *ex_what = nullptr;
    r_ext_ptr<mat> res;
    try {
        auto x = get<r_ext_ptr<mat>>(x_r);
        int k = get<int>(k_r);
        init_method im = parse_init_method(get<const char *>(method_r));
        switch (im) {
            case init_method::KMEANSPP:
                res.init(kmeanspp_init().init(*x, k));
                break;
            case init_method::RANDOM:
                res.init(random_init().init(*x, k));
                break;
            default:
                throw invalid_init_method("invalid init method");
        }
    } catch (exception &ex) {
        ex_what = ex.what();
    }

    if (ex_what)
        error(ex_what);

    try {
        SEXP r_res = put(*res);
        return r_res;
    } catch (exception &ex) {
        ex_what = ex.what();
        error(ex_what);
    }
}

R_CallMethodDef methods[] = {
        {"cec_r",              (DL_FUNC) &cec_r,              4},
        {"cec_split_r",        (DL_FUNC) &cec_split_r,        5},
        {"cec_init_centers_r", (DL_FUNC) &cec_init_centers_r, 3},
        {NULL, NULL,                                          0}
};

extern "C"
void R_init_CEC(DllInfo *dllInfo) {
    R_registerRoutines(dllInfo, NULL, methods, NULL, NULL);
    R_useDynamicSymbols(dllInfo, TRUE);
}
