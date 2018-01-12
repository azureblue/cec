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
                       int max_iter, int min_card)
                : x(x),
                  models(std::move(models)),
                  init(std::move(init)),
                  starts(starts),
                  max_iter(max_iter),
                  min_card(min_card) {}

        unique_ptr<clustering_results> operator()();

    private:
        const mat &x;
        vector<unique_ptr<model>> models;
        unique_ptr<initializer> init;
        int starts;
        int max_iter;
        int min_card;
        cec_starter starter;
        best_results_collector best;
    };

    class multi_starter {
    public:
        unique_ptr<clustering_results>
        start(const mat &x, vector<shared_ptr<model_spec>> models, init_spec init,
              int max_iter, int min_card, int starts, int threads);
    };
}

#endif
