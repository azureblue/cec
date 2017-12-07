#ifndef CEC_R_UTILS_H
#define    CEC_R_UTILS_H

#include "vec.h"
#include "exceptions.h"

extern "C" {
#include <Rinternals.h>
#include <Rdefines.h>
};

namespace cec {
    namespace r {

        template<typename T>
        T get(SEXP);

        template<>
        mat get<mat>(SEXP);

        template<>
        std::string get<std::string>(SEXP);

        template<>
        int get<int>(SEXP);

        template<>
        SEXP get<SEXP>(SEXP);

        template<>
        std::vector<int> get<std::vector<int>>(SEXP);

        template<>
        std::vector<double> get<std::vector<double>>(SEXP);

        SEXP put(const mat &ma);

        SEXP put(int);

        SEXP put(double);

        SEXP put(std::vector<int>);

        template<typename T>
        T get_named(SEXP list, const std::string &name) {
            int len = LENGTH(list);
            SEXP elementNames = GET_NAMES(list);
            for (int i = 0; i < len; i++) {
                if (name != CHAR(STRING_ELT(elementNames, i)))
                    continue;
                SEXP res = VECTOR_ELT(list, i);
                if (!res || isNull(res))
                    break;
                return cec::r::get<T>(res);
            }
            throw missing_parameter(name);
        };
        
        template<typename T>
        T get_n(SEXP vec, int idx){
            int len = LENGTH(vec);
            if (idx >= len)
                throw missing_parameter("out of range: " + idx);
            SEXP el = VECTOR_ELT(vec, idx);
            return cec::r::get<T>(el);
        }

        class r_wrapper {
            SEXP const sexp;

        public:
            r_wrapper(SEXP sexp) : sexp(sexp) {}

            int size() {
                return LENGTH(sexp);
            }

            template<typename T>
            T get(const std::string &name) {
                return get_named<T>(sexp, name);
            }

            template<typename T>
            T get() {
                return cec::r::get<T>(sexp);
            }

            template<typename T>
            T get(int idx) {
                return get_n<T>(sexp, idx);
            }
        };

    };
};
#endif    /* CEC_R_UTILS_H */

