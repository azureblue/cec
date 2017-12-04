#ifndef CEC_EXCEPTIONS_H
#define CEC_EXCEPTIONS_H

#include <exception>
#include "vec.h"
#include "cluster.h"

namespace cec {
    class invalid_covariance: public std::exception {
    public:
        explicit invalid_covariance(const cluster &cl, int cluster_number):
                cov(cl.covariance()),
                mean(cl.mean()),
                card(cl.card()),
                number(cluster_number){}
        const mat cov;
        const vec mean;
        const int card;
        const int number;
    };
}
#endif //CEC_EXCEPTIONS_H
