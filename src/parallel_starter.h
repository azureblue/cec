#ifndef CEC_PARALLEL_STARTER_H
#define CEC_PARALLEL_STARTER_H

#include <thread>
#include <future>
#include <utility>
#include "starter.h"
#include "multi_try_starter.h"
#include "split_starter.h"


namespace cec {
    template<class Task>
    class parallel_starter {
        using subtask = typename Task::subtask;
    public:
        parallel_starter(int threads, int starts)
                : threads(threads),
                  starts(starts) {}

        unique_ptr<clustering_results> start(Task task) {
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
                std::packaged_task<unique_ptr<clustering_results>()> task(
                        (starts_per_thread + (remaining-- > 0 ? 1 : 0)));
                cl_results.emplace_back(task.get_future());
                cl_tasks.emplace_back(std::move(task));
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
}

#endif //CEC_PARALLEL_STARTER_H
