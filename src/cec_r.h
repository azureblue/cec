#ifndef CEC_R_H
#define CEC_R_H

#include <Rinternals.h>
#include<R_ext/Rdynload.h>

extern "C" {
SEXP cec_r(SEXP x_r, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r);
SEXP cec_init_centers_r(SEXP x_r, SEXP k_r, SEXP method_r);
SEXP cec_split_r(SEXP x_r, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r,
                 SEXP split_param_r);
void R_init_CEC(DllInfo *dllInfo);
}

#endif /* CEC_R_H */
