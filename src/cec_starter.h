#include "starter.h"
#include "init.h"
#include "params.h"

#ifndef CEC_STARTER_H
#define CEC_STARTER_H

namespace cec {
    struct clustering_input {
        const mat &x;
        const vector<unique_ptr<model>> &models;

    public:
        clustering_input(const mat &x, const vector<unique_ptr<model>> &models)
                : x(x),
                  models(models) {}
    };

    class best_results_collector {
    public:
        void operator()(unique_ptr<clustering_results> &&cr) {
            if (!cr)
                return;
            if (!best || cr->energy < best->energy)
                best = std::move(cr);
        }

        unique_ptr<clustering_results> operator()() noexcept {
            return std::move(best);
        }

        void reset() {
            best.reset(nullptr);
        }

        explicit operator bool() {
            return !!best;
        }

    private:
        unique_ptr<clustering_results> best;
    };


    class clustering_starter {
    public:
        virtual ~clustering_starter() = default;
        virtual unique_ptr<clustering_results> start(const clustering_input &input_params) = 0;
    };

    class clustering_processor {
    public:
        virtual ~clustering_processor() = default;
        virtual unique_ptr<clustering_results>
        start(const unique_ptr<clustering_results> &cl_res, const clustering_input &input_param) = 0;
    };

    class cec_starter : public clustering_starter {
    public:
        struct parameters {
            const cec_parameters start_params;
            const centers_init_spec &init;
            int starts;

            parameters(cec_parameters start_params,
                       const centers_init_spec &init,
                       int starts = 1)
                    : start_params(start_params),
                      init(init),
                      starts(starts) {}
        };

        explicit cec_starter(const parameters &params)
                : starts(params.starts),
                  cec(params.start_params),
                  closest(),
                  init(params.init.create()) {}

        unique_ptr<clustering_results>
        start(const clustering_input &ip) override;

    private:
        const int starts;
        best_results_collector best;
        cross_entropy_clustering cec;
        closest_assignment closest;
        unique_ptr<centers_init> init;
    };
}

#endif
