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
        unique_ptr<clustering_results> start(Task &task, const mat &x,
                                             const vector<shared_ptr<model_spec>> &model_specs) {

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
                std::packaged_task<unique_ptr<clustering_results>(clustering_input &&)> pt(
                        task(starts_per_thread + (remaining-- > 0 ? 1 : 0)));
                cl_results.emplace_back(pt.get_future());
                cl_tasks.emplace_back(std::move(pt),
                                      clustering_input(x, model_spec::create_models(model_specs)));
            }

            best_results_collector best;

            best(my(clustering_input(x, model_spec::create_models(model_specs))));

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


    class mp_start_subtask {
    public:
        mp_start_subtask(mp_start_subtask &&) = default;

        mp_start_subtask(mp_start_subtask &) = delete;

        mp_start_subtask(unique_ptr<clustering_starter> c_starter,
                         vector<unique_ptr<clustering_processor>> &&c_procs, const int starts)
                : c_starter(std::move(c_starter)),
                  c_procs(std::move(c_procs)),
                  starts(starts) {};

        unique_ptr<clustering_results> operator()(clustering_input &&input_params) {
            best_results_collector best;
            for (int i = 0; i < starts; i++) {
                unique_ptr<clustering_results> res = c_starter->start(input_params);
                for (auto &&cp : c_procs)
                    res = cp->start(res, input_params);

                best(std::move(res));
            }
            return best();
        }

    private:
        unique_ptr<clustering_starter> c_starter;
        vector<unique_ptr<clustering_processor>> c_procs;
        const int starts;
    };

    class multiple_starts_task {
    public:
        using subtask = mp_start_subtask;

        explicit multiple_starts_task(cec_starter::parameters params)
                : params(params) {}

        subtask operator()(int starts) {
            unique_ptr<clustering_starter> starter = make_unique<cec_starter>(params);
            return mp_start_subtask(
                    std::move(starter), vector<unique_ptr<clustering_processor>>(), starts);
        }

    private:
        cec_starter::parameters params;
    };

    class start_and_split_task {
    public:
        using subtask = mp_start_subtask;

        explicit start_and_split_task(cec_starter::parameters init_cl_params,
                                      split_starter::parameters split_params
        )
                : init_cl_params(init_cl_params),
                  split_params(split_params) {}

        subtask operator()(int starts) {
            unique_ptr<clustering_starter> starter = make_unique<cec_starter>(init_cl_params);
            vector<unique_ptr<clustering_processor>> cl_procs(1);
            cl_procs[0] = make_unique<split_starter>(split_params);

            return mp_start_subtask(
                    std::move(starter),
                    std::move(cl_procs),
                    starts);
        }

    private:
        cec_starter::parameters init_cl_params;
        split_starter::parameters split_params;
    };
}

#endif //CEC_PARALLEL_STARTER_H
