#include "cec_starter.h"

namespace cec {
    unique_ptr<clustering_results> cec_starter::start(const clustering_input &ip) {
        const mat &x = ip.x;
        const vector<unique_ptr<model>> &models = ip.models;
        int k = models.size();
        best.reset();
        for (int i = 0; i < starts; i++)
            best(cec.start(x, closest.init(x, init->init(x, k)), models));
        return best();
    }
}
