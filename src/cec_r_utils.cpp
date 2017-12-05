#include "cec_r_utils.h"
#include "exceptions.h"

extern "C" {
#include <Rdefines.h>
};

template<typename T> std::string cec::r::get<std::string>(SEXP sexp) {
    if (isString(sexp))
        return std::string(CHAR(STRING_ELT(sexp, 0)));
    else
        throw invalid_parameter_type("string vector");
}
template<> cec::mat cec::r::get<cec::mat>(SEXP R_ma) {

    int m = Rf_nrows(R_ma);
    int n = Rf_ncols(R_ma);
    mat ma(m, n);
    double *m_data = ma.data();
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            m_data[i * m + j] = REAL(R_ma)[j * m + i];
    return ma;
}

SEXP cec::r::create_R_matrix(const cec::mat &ma) {
    int m = ma.m;
    int n = ma.n;
    const double *m_data = ma.data();

    SEXP R_ma;

    PROTECT(R_ma = allocMatrix(REALSXP, m, n));

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            REAL(R_ma)[j * m + i] = m_data[i * m + j];

    UNPROTECT(1);

    return R_ma;
}

template <typename T>
static T cec::r::get_named_element(SEXP list, std::string &name) {
    int len = LENGTH(list);
    SEXP elementNames = GET_NAMES(list);
    for (int i = 0; i < len; i++)
        if (name == CHAR(STRING_ELT(elementNames, i))) {
            SEXP res = VECTOR_ELT(list, i);
            if (!res || isNull(res))
                break;

            return get(res);
        }
    throw missing_parameter(name);
};
