#ifndef CEC_R_H
#define CEC_R_H

extern "C" {
#include <Rinternals.h>
SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r);
}

#endif /* CEC_R_H */
