#ifndef CEC_EXCEPTIONS_H
#define CEC_EXCEPTIONS_H

#include <exception>
#include <utility>
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

        const char *what() const noexcept override {
            return "invalid covariance: probably not positive definite";
        }
    };

    class all_clusters_removed: public std::exception {
    public:
        const char *what() const noexcept override {
            return "all clusters have been removed";
        }
    };

    class not_implemented: public std::exception {
    public:
        const std::string desc;

        explicit not_implemented(std::string desc) : desc(std::move(desc)) {}

        const char *what() const noexcept override {
            return ("not implemented: " + desc).c_str();
        }
    };

    class invalid_init_method: public std::exception {
    public:
        const std::string method;

        explicit invalid_init_method(std::string method): method(std::move(method)) {}

        const char *what() const noexcept override {
            return ("invalid center initialization method: " + method).c_str();
        }
    };

    class invalid_model_name: public std::exception {
    public:
        const std::string name;

        explicit invalid_model_name(std::string name): name(std::move(name)) {}

        const char *what() const noexcept override {
            return ("invalid model name: " + name).c_str();
        }
    };

    class missing_parameter: public std::exception {
    public:
        const std::string name;

        explicit missing_parameter(std::string name): name(std::move(name)) {}

        const char *what() const noexcept override {
            return ("missing parameter: " + name).c_str();
        }
    };

    class invalid_parameter_type: public std::exception {
    public:
        const std::string expected;

        explicit invalid_parameter_type(std::string name): expected(std::move(name)) {}

        const char *what() const noexcept override {
            return ("invalid parameter type, expected: " + expected).c_str();
        }
    };
}
#endif //CEC_EXCEPTIONS_H
