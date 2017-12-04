#ifndef CEC_R_UTILS_H
#define	CEC_R_UTILS_H

#include "vec.h"

extern "C" {
#include <Rinternals.h>
};

namespace cec::r {
    template<typename T> T get(SEXP sexp);
    SEXP create_R_matrix(const mat &);
    template <typename T>
    T get_named_element(SEXP list, std::string name);
}
#endif	/* CEC_R_UTILS_H */

