#ifndef CEC_PARALLEL_STARTER_H
#define CEC_PARALLEL_STARTER_H

#include <thread>
#include <future>
#include <utility>
#include "starter.h"
#include "cec_starter.h"
#include "split_starter.h"


namespace cec {
    class parallel_starter {
    public:
        parallel_starter(int threads, int starts)
                : threads(threads),
                  starts(starts) {}

        template<class Task>
        unique_ptr<clustering_results> start(Task task) {
            using subtask = typename Task::subtask;

            if (threads == 0)
                threads = std::thread::hardware_concurrency();
            if (threads == 0)
                threads = default_threads_number;
            threads = std::min(starts, threads);
            int starts_per_thread = starts / threads;
            int remaining = starts - (starts_per_thread * threads);

            vector<std::thread> cl_tasks;
            vector<std::future<unique_ptr<clustering_results>>> cl_results;

            subtask my = task(starts_per_thread + (remaining-- > 0 ? 1 : 0));

            for (int th = 1; th < threads; th++) {
                std::packaged_task<unique_ptr<clustering_results>()> pt(
                        task(starts_per_thread + (remaining-- > 0 ? 1 : 0)));
                cl_results.emplace_back(pt.get_future());
                cl_tasks.emplace_back(std::move(pt));
            }

            best_results_collector best;

            best(my());

            for (auto &&cl_task : cl_tasks)
                cl_task.join();

            for (auto &&result : cl_results)
                best(result.get());

            return best();
        }

    private:
        int threads;
        int starts;
        static const int default_threads_number = 4;
    };

    class multiple_starts_task {
    public:
        class mp_start_subtask {
        public:
            mp_start_subtask(mp_start_subtask &&) = default;
            mp_start_subtask(mp_start_subtask &) = delete;

            mp_start_subtask(cec_starter::parameters params, const mat &x, const vector <shared_ptr<model_spec>> &models)
                    : ms(params),
                      x(x),
                      models(models) {}

            unique_ptr<clustering_results> operator()() {
                return ms.start(x, models);
            }

        private:
            cec_starter ms;
            const mat &x;
            const vector<shared_ptr<model_spec>> &models;
        };

        using subtask = mp_start_subtask;

        multiple_starts_task(const mat &x, const vector <shared_ptr<model_spec>> &models, cec_parameters start_params,
                                  const centers_init_spec &init_spec)
                :  x(x),
                  models(models),
                  init_spec(init_spec),
                  cec_params(start_params) {}


        subtask operator()(int starts) {
            return mp_start_subtask(cec_starter::parameters(cec_params, init_spec, starts), x, models);
        }

    private:
        const mat &x;
        const vector<shared_ptr<model_spec>> &models;
        const centers_init_spec &init_spec;
        cec_parameters cec_params;
    };
}

#endif //CEC_PARALLEL_STARTER_H
