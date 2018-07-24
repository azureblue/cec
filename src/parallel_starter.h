#ifndef CEC_PARALLEL_STARTER_H
#define CEC_PARALLEL_STARTER_H

#include <thread>
#include <future>
#include <utility>
#include "starter.h"
#include "cec_starter.h"
#include "split_starter.h"
#include "exceptions.h"

namespace cec {

    class unique_models_input {
    public:
        unique_models_input(const mat &x, vector<unique_ptr<model>> &&models)
                : x(x),
                  models(std::move(models)) {}

        unique_models_input(unique_models_input &) = delete;

        unique_models_input(unique_models_input &&) = default;

        clustering_input get() {
            return clustering_input(x, models);
        }

    private:
        const mat &x;
        vector<unique_ptr<model>> models;
    };

    class parallel_starter {
    public:
        parallel_starter(int max_threads, int starts)
                : starts(starts) {
            if (max_threads == 0)
                max_threads = std::thread::hardware_concurrency();
            if (max_threads == 0)
                max_threads = default_threads_number;

            parallel_starter::max_threads = std::min(starts, max_threads);
        }

        template<class Task>
        unique_ptr<clustering_results> start(Task &task) {
            using subtask = typename Task::subtask;
            int starts_per_thread = starts / max_threads;
            int remaining = starts - (starts_per_thread * max_threads);
            best_results_collector best;

            vector<std::thread> cl_tasks;
            vector<std::future<unique_ptr<clustering_results>>> cl_results;

            subtask my = task(starts_per_thread + (remaining-- > 0 ? 1 : 0));
            try {
                for (int th = 1; th < max_threads; th++) {
                    std::packaged_task<unique_ptr<clustering_results>()>
                            pt_subtask(task(starts_per_thread + (remaining-- > 0 ? 1 : 0)));

                    cl_results.emplace_back(pt_subtask.get_future());
                    cl_tasks.emplace_back(std::move(pt_subtask));
                }

                best(my());

                for (auto &&result : cl_results)
                    try {
                        best(result.get());
                    } catch (clustering_exception &ex) {
                        //ignore for now...
                    }

            } catch (std::exception &ex) {
                join_all_threads(cl_tasks);
                throw;
            }

            join_all_threads(cl_tasks);

            return best();
        }

    private:
        int max_threads;
        int starts;
        static const int default_threads_number = 4;

        static void join_all_threads(vector<std::thread> &threads) {
            for (auto &&th : threads)
                th.join();
        }
    };

    class mp_start_subtask {
    public:
        mp_start_subtask(mp_start_subtask &&) = default;

        mp_start_subtask(mp_start_subtask &) = delete;

        mp_start_subtask(unique_ptr<clustering_starter> c_starter,
                         vector<unique_ptr<clustering_processor>> c_procs, unique_models_input &&uniqe_m_input,
                         const int starts)
                : c_starter(std::move(c_starter)),
                  c_procs(std::move(c_procs)),
                  uniqe_m_input(std::move(uniqe_m_input)),
                  starts(starts) {};

        unique_ptr<clustering_results> operator()() {
            best_results_collector best;
            for (int i = 0; i < starts; i++) {
                unique_ptr<clustering_results> res = c_starter->start(uniqe_m_input.get());
                for (auto &&cp : c_procs)
                    res = cp->start(res, uniqe_m_input.get());

                best(std::move(res));
            }
            return best();
        }

    private:
        unique_ptr<clustering_starter> c_starter;
        vector<unique_ptr<clustering_processor>> c_procs;
        unique_models_input uniqe_m_input;
        const int starts;
    };

    class multiple_starts_task {
    public:
        using subtask = mp_start_subtask;

        explicit multiple_starts_task(cec_starter::parameters params, const mat &x,
                                      const vector<shared_ptr<model_spec>> &model_specs)
                : params(params),
                  x(x),
                  model_specs(model_specs) {}

        subtask operator()(int starts) const {
            unique_ptr<clustering_starter> starter = make_unique<cec_starter>(params);
            return mp_start_subtask(std::move(starter),
                                    vector<unique_ptr<clustering_processor>>(),
                                    unique_models_input(x, model_spec::create_models(model_specs)), starts);
        }

    private:
        const cec_starter::parameters params;
        const mat &x;
        const vector<shared_ptr<model_spec>> &model_specs;
    };

    class start_and_split_task {
    public:
        using subtask = mp_start_subtask;

        explicit start_and_split_task(cec_starter::parameters init_cl_params,
                                      split_starter::parameters split_params,
                                      const mat &x,
                                      const vector<shared_ptr<model_spec>> &model_specs
        )
                : init_cl_params(init_cl_params),
                  split_params(split_params),
                  x(x),
                  model_specs(model_specs) {}

        subtask operator()(int starts) {
            unique_ptr<clustering_starter> starter = make_unique<cec_starter>(init_cl_params);
            vector<unique_ptr<clustering_processor>> cl_procs(1);
            cl_procs[0] = make_unique<split_starter>(split_params);

            return mp_start_subtask(
                    std::move(starter),
                    std::move(cl_procs),
                    unique_models_input(x, model_spec::create_models(model_specs)),
                    starts);
        }

    private:
        cec_starter::parameters init_cl_params;
        split_starter::parameters split_params;
        const mat &x;
        const vector<shared_ptr<model_spec>> &model_specs;
    };
}

#endif //CEC_PARALLEL_STARTER_H
