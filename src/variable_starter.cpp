#include "variable_starter.h"

namespace cec {
    unique_ptr<clustering_results>
    variable_starter::start(const mat &x, vector<shared_ptr<model_spec>> m_specs) {
        best_results_collector best;
        for (auto &&k : centers_number) {
            vector<shared_ptr<model_spec>> models_specs_subset(m_specs.begin(), m_specs.begin() + k);
            try {
                best(cl_starter(x, models_specs_subset));
            } catch (clustering_exception &ce) {}
        }
        return best();
    }
}
