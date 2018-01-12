#include "variable_starter.h"

std::unique_ptr<cec::clustering_results>
cec::variable_starter::start(const mat &x, model_specs m_specs, initializer_spec init_spec,
                             const variable_starter_params &params) {
    best_results_collector best;
    for (auto &&k : params.centers_number)
        best(ms.start(x, m_specs, init_spec, params.msp));
    return best();
}
