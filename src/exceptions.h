#ifndef CEC_EXCEPTIONS_H
#define CEC_EXCEPTIONS_H

#include <exception>
#include <utility>
#include "common.h"
#include "vec.h"
#include "cluster.h"

namespace cec {
    class cec_exception : public std::exception {
    public:
        const string &info() {
            return info_msg;
        }

    protected:
        string info_msg;

        explicit cec_exception(string info)
                : info_msg(std::move(info)) {}

    public:
        const char *what() const noexcept override {
            return "cec exception";
        }
    };

    class clustering_exception : public cec_exception {
    public:
        explicit clustering_exception(string info)
                : cec_exception(std::move(info)) {}

        const char *what() const noexcept override {
            return "clustering failed";
        }
    };

    class invalid_covariance : public clustering_exception {
    public:
        explicit invalid_covariance(mat cov)
                : clustering_exception("invalid covariance: probably not positive definite"),
                  cov(std::move(cov)) {}

        const mat &covariance() const {
            return cov;
        }

    private:
        mat cov;
    };


    class all_clusters_removed : public clustering_exception {
    public:
        all_clusters_removed()
                : clustering_exception("all clusters have been removed") {}
    };

    class not_implemented : public cec_exception {
    public:
        explicit not_implemented(const string &info)
                : cec_exception(info) {}

        const char *what() const noexcept override {
            return "not implemented";
        }
    };

    class invalid_init_method : public cec_exception {
    public:
        explicit invalid_init_method(string info)
                : cec_exception(std::move(info)) {}

        const char *what() const noexcept override {
            return "invalid center method";
        }
    };

    class invalid_model_name : public cec_exception {
    public:
        explicit invalid_model_name(string name)
                : cec_exception(std::move(name)) {}

        const char *what() const noexcept override {
            return "invalid model name";
        }
    };

    class missing_parameter : public cec_exception {
    public:
        explicit missing_parameter(std::string name)
                : cec_exception(std::move(name)) {}

        const char *what() const noexcept override {
            return "missing parameter";
        }
    };

    class invalid_parameter_type : public cec_exception {
    public:
        explicit invalid_parameter_type(string name)
                : cec_exception(std::move(name)) {}

        const char *what() const noexcept override {
            return "invalid parameter type";
        }
    };

    class invalid_model_parameter : public cec_exception {
    public:
        explicit invalid_model_parameter(string desc)
                : cec_exception(std::move(desc)) {}

        const char *what() const noexcept override {
            return "invalid model parameter";
        }
    };
}
#endif //CEC_EXCEPTIONS_H
