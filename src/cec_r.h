#ifndef CEC_R_H
#define CEC_R_H

#include <Rinternals.h>
extern "C" {
SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r);
SEXP cec_init_centers_r(SEXP x_r, SEXP k_r, SEXP method_r);
}

#endif /* CEC_R_H */
