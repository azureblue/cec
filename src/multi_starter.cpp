#include "multi_starter.h"
#include <thread>
#include <future>

namespace cec {
    unique_ptr<clustering_results> clustering_task::operator()() {
        best_results_collector best;
        closest_assignment ca;
        int k = models.size();
        for (int i = 0; i < starts; i++)
            best(starter.start(x, ca.init(x, init->init(x, k)), models, start_params));
        return best();
    }

    unique_ptr<clustering_results>
    multi_starter::start(const mat &x,
                         vector<shared_ptr<model_spec>> models,
                         const centers_init_spec &init,
                         const multi_starter_params &param) {
        int threads = param.threads;
        int starts = param.starts;
        if (threads == 0)
            threads = std::thread::hardware_concurrency();
        if (threads == 0)
            threads = default_threads_number;
        threads = std::min(starts, threads);
        int starts_per_thread = starts / threads;
        int remaining = starts - (starts_per_thread * threads);

        vector<std::thread> cl_tasks;
        vector<std::future<unique_ptr<clustering_results>>> cl_results;

        clustering_task my_task(x, model_spec::create_models(models), init.create(),
                                starts_per_thread + (remaining-- > 0 ? 1 : 0), param.start_params);

        for (int th = 1; th < threads; th++) {
            std::packaged_task<unique_ptr<clustering_results>()> task(
                    clustering_task(x, model_spec::create_models(models), init.create(),
                                    starts_per_thread + (remaining-- > 0 ? 1 : 0),
                                    param.start_params));
            cl_results.emplace_back(task.get_future());
            cl_tasks.emplace_back(std::move(task));
        }

        best_results_collector best;

        best(my_task());

        for (auto &&cl_task : cl_tasks)
            cl_task.join();

        for (auto &&result : cl_results)
            best(result.get());

        return best();
    }
}
