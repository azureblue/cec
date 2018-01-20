#include "cec_starter.h"

namespace cec {
    unique_ptr<clustering_results>
    cec_starter::start(const mat &x, const vector<shared_ptr<model_spec>> &m_specs) {
        int starts = params.starts;
        best_results_collector best;
        int k = m_specs.size();
        vector<unique_ptr<model>> models = model_spec::create_models(m_specs);
        for (int i = 0; i < starts; i++)
            best(cec.start(x, closest.init(x, init->init(x, k)), models));

        return best();
    }
}
