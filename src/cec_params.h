#ifndef CEC_PARAMS_H
#define CEC_PARAMS_H

#include <stdbool.h>
#include "alloc.h"
#include "matrix.h"

enum centers_init_method {
    NONE = 0,
    KMEANSPP = 1,
    RANDOM = 2
};

enum density_family {
    ALL, COVARIANCE, DIAGONAL, EIGENVALUES, FIXED_R, SPHERICAL
};

bool parse_init_method(const char * method, enum centers_init_method * result);
bool cec_parse_type(const char * type, enum density_family * result);

struct cec_centers_param {
    enum centers_init_method init_m;
    cec_mat * centers_mat;
    const vec_i * var_centers;
};

struct cec_control_param {
    int starts;
    int max_iterations;
    int min_card;
    int threads;
};

struct cec_model_r_params {
    double r;
};

struct cec_model_eigenvalues_params {
    const vec_d * given_eigenvalues;
};

struct cec_model_covariances_params {
    const cec_mat * cov;
    const cec_mat * cov_inv;
};

struct cec_model_specification {
    enum density_family type;
    int n;
    union {
        struct cec_model_r_params r_params;
        struct cec_model_eigenvalues_params eigenvalues_params;
        struct cec_model_covariances_params covariances_params;
    };
};

struct cec_models_param {
    int len;
    struct cec_model_specification *model_specs;
};

typedef struct cec_centers_param cec_centers_par;
typedef struct cec_control_param cec_control_par;
typedef struct cec_models_param cec_models_par;
typedef struct cec_model_specification cec_model_spec;

#endif //CEC_PARAMS_H
