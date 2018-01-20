#ifndef CEC_SPLIT_STARTER_H
#define CEC_SPLIT_STARTER_H

#include "starter.h"
#include "cec_starter.h"

namespace cec {

    class split_starter {
        struct parameters {
            cec_parameters start_params;
            const centers_init_spec &init;
            int split_tries;
            int max_k;
            int max_depth;
        };

    public:
        split_starter(const parameters &params)
                : params(params),
                  splitter(cec_starter::parameters(params.start_params, params.init,
                                                        params.split_tries)),
                  cec(params.start_params) {}

        unique_ptr<clustering_results>
        start(const unique_ptr<clustering_results> &cl_res, const mat &x_mat,
              const shared_ptr<model_spec> &model_sp);

    private:
        parameters params;
        cec_starter splitter;
        cross_entropy_clustering cec;

        unique_ptr<clustering_results>
        try_split_cluster(const mat &x_mat, shared_ptr<model_spec> m_spec);


    };
}

#endif //CEC_SPLIT_STARTER_H
