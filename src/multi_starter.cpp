#include "multi_starter.h"

std::unique_ptr<cec::clustering_results> cec::clustering_task::operator()() {
    int k = models.size();
    for (int i = 0; i < starts; i++) {
        best(starter.start(x, init->init(x, k), models, start_params));
    }
    return best();
}

std::unique_ptr<cec::clustering_results>
cec::multi_starter::start(const cec::mat &x, std::vector<std::shared_ptr<cec::model_spec>> models,
                          const initializer_spec &init,
                          const multi_starter_params &param) {
    int threads = param.threads;
    int starts = param.starts;
    if (threads == 0)
        threads = std::thread::hardware_concurrency();
    if (threads == 0)
        threads = 4;
    threads = std::min(starts, threads);
    int starts_per_thread = starts / threads;
    int remaining = starts - (starts_per_thread * threads);

    vector<std::thread> cl_tasks;
    vector<std::future<unique_ptr<clustering_results>>> cl_results;

    clustering_task my_task(x, model_spec::create_models(models), init.create(),
                            starts_per_thread + (remaining-- > 0 ? 1 : 0), param.start_params
    );

    for (int th = 1; th < threads; th++) {
        std::packaged_task<unique_ptr<clustering_results>()> task(
                clustering_task(x, model_spec::create_models(models), init.create(),
                                starts_per_thread + (remaining-- > 0 ? 1 : 0), param.start_params
                ));
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
