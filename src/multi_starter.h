#include "starter.h"
#include "init.h"
#include "params.h"

#ifndef CEC_MULTI_STARTER_H
#define CEC_MULTI_STARTER_H

namespace cec {
    class best_results_collector {
        unique_ptr<clustering_results> best;
    };
    class clusterig_thread {
    public:
        void operator()() {
            int k = models.size();
            for (int i = 0; i < starts; i++) {
                const vector<int> &asgn = a_init->init(x, c_init->init(x, k));
                const clustering_results &res = starter.start(x, asgn, models, max_iter, min_card);
            }
        }

    private:
        const mat &x;
        const vector<unique_ptr<model>> models;
        const unique_ptr<centers_init> c_init;
        const unique_ptr<assignment_init> a_init;
        const int starts;
        const int max_iter;
        const int min_card;
        cec_starter starter;

    };
    class multi_starter {
    public:
        clustering_results start(const mat &x, const vector<shared_ptr<model_spec>> &models,
                                 int max_iter, int min_card, int threads) {

        }
    };
}

#endif
