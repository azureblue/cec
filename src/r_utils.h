#ifndef CEC_R_UTILS_H
#define CEC_R_UTILS_H

#include "vec.h"
#include "exceptions.h"
#include "r_ext_ptr.h"

#include <Rinternals.h>
#include <Rdefines.h>

namespace cec {

    namespace r {
        template<typename T>
        inline T get(SEXP sexp);

        template<>
        inline const char* get(SEXP sexp) {
            if (!isString(sexp))
                throw invalid_parameter_type("string vector");
            return CHAR(STRING_ELT(sexp, 0));
        }

        template<>
        inline int get(SEXP sexp) {
            if (TYPEOF(sexp) != INTSXP || LENGTH(sexp) != 1)
                throw invalid_parameter_type("single integer");
            return INTEGER(sexp)[0];
        }

        template<>
        inline double get(SEXP sexp) {
            if (TYPEOF(sexp) != REALSXP || LENGTH(sexp) != 1)
                throw invalid_parameter_type("single real");
            return REAL(sexp)[0];
        }

        template<>
        inline r_ext_ptr<mat> get(SEXP sexp) {
            if (!isMatrix(sexp))
                throw invalid_parameter_type("matrix");
            int m = Rf_nrows(sexp);
            int n = Rf_ncols(sexp);
            double *r_ma_data = REAL(sexp);
            auto ma = make_r_ext<mat>(m, n);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    (*ma)[i][j] = r_ma_data[j * m + i];

            return ma;
        }

        template<>
        inline r_ext_ptr<vector<double>> get(SEXP sexp) {
            if (TYPEOF(sexp) != REALSXP)
                throw invalid_parameter_type("real vector");
            return make_r_ext<vector<double>>(REAL(sexp), REAL(sexp) + LENGTH(sexp));
        }

        template<>
        inline r_ext_ptr<vector<int>> get(SEXP sexp) {
            if (TYPEOF(sexp) != INTSXP)
                throw invalid_parameter_type("integer vector");
            return make_r_ext<vector<int>>(INTEGER(sexp), INTEGER(sexp) + LENGTH(sexp));
        }

        inline SEXP put(const mat &ma) {
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

        inline SEXP put(int val) {
            SEXP ve;
            PROTECT(ve = allocVector(INTSXP, 1));
            INTEGER(ve)[0] = val;
            UNPROTECT(1);

            return ve;
        }

        inline SEXP put(double val) {
            SEXP ve;
            PROTECT(ve = allocVector(REALSXP, 1));
            REAL(ve)[0] = val;
            UNPROTECT(1);

            return ve;
        }

        inline SEXP put(vector<int> val) {
            SEXP ve;
            PROTECT(ve = allocVector(INTSXP, val.size()));
            std::copy(val.begin(), val.end(), INTEGER(ve));
            UNPROTECT(1);

            return ve;
        }

        class r_wrapper {
        public:
            explicit r_wrapper(SEXP sexp)
                    : sexp(sexp) {}

            int size() {
                return LENGTH(sexp);
            }

            r_wrapper operator[](const char *name) {
                return r_wrapper(get_named(sexp, name));
            }

            r_wrapper operator[](const int idx) {
                return r_wrapper(get_n(sexp, idx));
            }

            template<typename T>
            T get() {
                return cec::r::get<T>(sexp);
            }

        private:
            SEXP sexp;

            SEXP get_named(SEXP list, const char *name) {
                SEXP elementNames = GET_NAMES(list);
                if (!isString(elementNames))
                    throw invalid_parameter_type("named elements");
                int len = LENGTH(elementNames);
                for (int i = 0; i < len; i++) {
                    if (strcmp(name, CHAR(STRING_ELT(elementNames, i))))
                        continue;
                    SEXP res = VECTOR_ELT(list, i);
                    if (!res || isNull(res))
                        break;
                    return res;
                }
                throw missing_parameter(name);
            }

            SEXP get_n(SEXP vec, int idx) {
                int len = LENGTH(vec);
                if (idx >= len)
                    throw missing_parameter("out of range: " + idx);
                SEXP el = VECTOR_ELT(vec, idx);
                return el;
            }
        };
    }
}
#endif    /* CEC_R_UTILS_H */

