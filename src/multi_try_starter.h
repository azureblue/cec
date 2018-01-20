#include "starter.h"
#include "init.h"
#include "params.h"

#ifndef CEC_MULTI_STARTER_H
#define CEC_MULTI_STARTER_H

namespace cec {
    class best_results_collector {
    public:
        void operator()(unique_ptr<clustering_results> &&cr) {
            if (!cr)
                return;
            if (!best || cr->energy < best->energy)
                best = std::move(cr);
        }

        unique_ptr<clustering_results> operator()() {
            return std::move(best);
        }

        explicit operator bool() {
            return !!best;
        }

    private:
        unique_ptr<clustering_results> best;
    };


    class clustering_starter {
    public:
        virtual unique_ptr<clustering_results> start(const mat &x, const vector<shared_ptr<model_spec>> &m_specs) = 0;
    };

    class multi_try_starter: public clustering_starter {
    public:
        struct parameters {
            starter_params start_params;
            const centers_init_spec &init;
            int starts;

            parameters(const starter_params &start_params, const centers_init_spec &init,
                       int starts)
                    : start_params(start_params),
                      init(init),
                      starts(starts) {}
        };

        explicit multi_try_starter(const parameters &params)
                : params(params),
                  cec(params.start_params),
                  closest(),
                  init(params.init.create()) {}

        unique_ptr<clustering_results>
        start(const mat &x, const vector<shared_ptr<model_spec>> &models);

    private:
        parameters params;
        cec_starter cec;
        closest_assignment closest;
        unique_ptr<centers_init> init;
    };
}

#endif
