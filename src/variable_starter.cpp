#include "variable_starter.h"

namespace cec {
    unique_ptr<clustering_results>
    variable_starter::start(const mat &x, model_specs m_specs, initializer_spec init_spec,
                            const variable_starter_params &params) {
        best_results_collector best;
        for (auto &&k : params.centers_number) {
            model_specs models_specs_subset(m_specs.begin(), m_specs.begin() + k);
            best(ms.start(x, models_specs_subset, init_spec, params.ms_params));
        }
        return best();
    }
}