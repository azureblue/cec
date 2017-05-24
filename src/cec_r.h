#ifndef CEC_R_H
#define	CEC_R_H

#include <Rinternals.h>
#include "centers_init.h"

/*
 * Entry point to CEC - called from R.
 */

struct cec_centers_param {
    init_method init_m;
    cec_mat * centers_mat;
    const int * var_centers;
    int var_centers_len;
};

struct cec_control_param {
    int starts;
    int max_iterations;
    int min_card;
};

typedef struct cec_centers_param cec_centers;
typedef struct cec_control_param cec_control;

SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP type, SEXP params);

#endif	/* CEC_R_H */
