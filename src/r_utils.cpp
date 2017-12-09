#include "r_utils.h"
#include "exceptions.h"

template<>
std::string cec::r::get<std::string>(SEXP sexp) {
    if (!isString(sexp))
        throw invalid_parameter_type("string vector");
    return std::string(CHAR(STRING_ELT(sexp, 0)));
}

template<>
SEXP cec::r::get<SEXP>(SEXP sexp) {
    return sexp;
}

template<>
int cec::r::get<int>(SEXP sexp) {
    if (TYPEOF(sexp) != INTSXP)
        throw invalid_parameter_type("integer");
    return asInteger(sexp);
}

template<>
cec::mat cec::r::get<cec::mat>(SEXP sexp) {
    if (!isMatrix(sexp))
        throw invalid_parameter_type("matrix");
    int m = Rf_nrows(sexp);
    int n = Rf_ncols(sexp);
    mat ma(m, n);
    double *m_data = ma.data();
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            m_data[i * n + j] = REAL(sexp)[j * m + i];
    return ma;
}

template<>
std::vector<int> cec::r::get<std::vector<int>>(SEXP sexp) {
    if (TYPEOF(sexp) != INTSXP)
        throw invalid_parameter_type("integer vector");
    int len = LENGTH(sexp);
    std::vector<int> vec_i(len);
    int *r_data = INTEGER(sexp);
    std::copy(&r_data[0], &r_data[len], vec_i.begin());
    return vec_i;
}

template<>
std::vector<double> cec::r::get<std::vector<double>>(SEXP sexp) {
    if (TYPEOF(sexp) != REALSXP)
        throw invalid_parameter_type("real vector");
    int len = LENGTH(sexp);
    std::vector<double> vec_d(len);
    double *r_data = REAL(sexp);
    std::copy(&r_data[0], &r_data[len], vec_d.begin());
    return vec_d;
}


SEXP cec::r::put(const cec::mat &ma) {
    int m = ma.m;
    int n = ma.n;
    const double *m_data = ma.data();

    SEXP r_ma;

    PROTECT(r_ma = allocMatrix(REALSXP, m, n));

    double *r_data = REAL(r_ma);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            r_data[j * m + i] = m_data[i * n + j];

    UNPROTECT(1);

    return r_ma;
}

SEXP cec::r::put(int val) {
    SEXP ve;
    PROTECT(ve = allocVector(INTSXP, 1));
    INTEGER(ve)[0] = val;
    UNPROTECT(1);

    return ve;
}

SEXP cec::r::put(double val) {
    SEXP ve;
    PROTECT(ve = allocVector(REALSXP, 1));
    REAL(ve)[0] = val;
    UNPROTECT(1);

    return ve;
}

SEXP cec::r::put(std::vector<int> val) {
    SEXP ve;
    PROTECT(ve = allocVector(INTSXP, val.size()));
    std::copy(val.begin(), val.end(), INTEGER(ve));
    UNPROTECT(1);

    return ve;
}
