#include "variable_starter.h"

namespace cec {
    unique_ptr<clustering_results>
    variable_starter::start(const mat &x, vector<shared_ptr<model_spec>> m_specs) {
        best_results_collector best;
        for (auto &&k : centers_number) {
            vector<shared_ptr<model_spec>> models_specs_subset(m_specs.begin(), m_specs.begin() + k);
            best(cs->start(clustering_input(x, model_spec::create_models(models_specs_subset))));
        }
        return best();
    }
}