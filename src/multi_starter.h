#include "starter.h"
#include "init.h"
#include "params.h"
#include <thread>
#include <future>

#ifndef CEC_MULTI_STARTER_H
#define CEC_MULTI_STARTER_H

namespace cec {
    class best_results_collector {
    public:
        void operator()(unique_ptr<clustering_results> &&cr) {
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

    class clustering_task {
    public:
        clustering_task(const mat &x, vector<unique_ptr<model>> &&models,
                       unique_ptr<initializer> &&init, int starts,
                       const starter_params &start_params)
                : x(x),
                  models(std::move(models)),
                  init(std::move(init)),
                  starts(starts),
                  start_params(start_params)
        {}

        unique_ptr<clustering_results> operator()();

    private:
        const mat &x;
        vector<unique_ptr<model>> models;
        unique_ptr<initializer> init;
        int starts;
        const starter_params &start_params;
        cec_starter starter;
        best_results_collector best;
    };

    struct multi_starter_params {
        starter_params start_params;
        int starts;
        int threads;

        multi_starter_params(const starter_params &start_params, int starts, int threads)
                : start_params(start_params),
                  starts(starts),
                  threads(threads) {}
    };

    class multi_starter {
    public:
        unique_ptr<clustering_results>
        start(const mat &x, vector<shared_ptr<model_spec>> models, const initializer_spec &init,
              const multi_starter_params &params);
    };
}

#endif
