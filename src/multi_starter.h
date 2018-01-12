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

        operator bool() {
            return !!best;
        }

    private:
        unique_ptr<clustering_results> best;
    };

    class clusterig_task {
    public:
        clusterig_task(const mat &x, vector<unique_ptr<model>> &&models, unique_ptr<initializer> &&init, int starts,
                       int max_iter, int min_card)
                : x(x),
                  models(std::move(models)),
                  init(std::move(init)),
                  starts(starts),
                  max_iter(max_iter),
                  min_card(min_card) {}

        unique_ptr<clustering_results> operator()() {
            int k = models.size();
            for (int i = 0; i < starts; i++) {
                br(starter.start(x, init->init(x, k), models, max_iter, min_card));
            }
            return br();
        }

    private:
        const mat &x;
        vector<unique_ptr<model>> models;
        unique_ptr<initializer> init;
        int starts;
        int max_iter;
        int min_card;
        cec_starter starter;
        best_results_collector br;

    };

    class multi_starter {
    public:
        unique_ptr<clustering_results> start(const mat &x, vector<shared_ptr<model_spec>> models, init_spec init,
                                             int max_iter, int min_card, int starts, int threads) {
            if (threads == 0)
                threads = std::thread::hardware_concurrency();
            if (threads == 0)
                threads = 4;
            threads = std::min(starts, threads);
            int starts_per_thread = starts / threads;
            int remaining = starts - (starts_per_thread * threads);

            vector<std::thread> cl_tasks;
            vector<std::future<unique_ptr<clustering_results>>> cl_results;

            clusterig_task my_task(x, model_spec::create_models(models), init.create(),
                           starts_per_thread + (remaining-- > 0 ? 1 : 0), max_iter, min_card
            );

            for (int th = 1; th < threads; th++) {
                std::packaged_task<unique_ptr<clustering_results>()> task(
                        clusterig_task(x, model_spec::create_models(models), init.create(),
                                       starts_per_thread + (remaining-- > 0 ? 1 : 0), max_iter, min_card
                        ));
                cl_results.emplace_back(task.get_future());
                cl_tasks.emplace_back(std::move(task));
            }

            best_results_collector br;

            br(my_task());

            for (auto &&cl_task : cl_tasks)
                cl_task.join();
            for (auto &&result : cl_results)
                br(result.get());

            return br();

        }
    };
}

#endif
