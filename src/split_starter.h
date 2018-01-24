#ifndef CEC_SPLIT_STARTER_H
#define CEC_SPLIT_STARTER_H

#include "starter.h"
#include "cec_starter.h"

namespace cec {

    class split_starter : public clustering_processor {
    public:
        struct parameters {
            cec_parameters start_params;
            const model_spec &model_sp;
            const centers_init_spec &init;
            int split_tries;
            int max_k;
            int max_depth;

            parameters(const cec_parameters &start_params, const model_spec &model_sp,
                       const centers_init_spec &init, int split_tries, int max_k, int max_depth)
                    : start_params(start_params),
                      model_sp(model_sp),
                      init(init),
                      split_tries(split_tries),
                      max_k(max_k),
                      max_depth(max_depth) {}
        };

        explicit split_starter(const parameters &params)
                : splitter(cec_starter::parameters(params.start_params, params.init,
                                                   params.split_tries)),
                  cec(params.start_params),
                  m_spec(params.model_sp),
                  max_k(params.max_k),
                  max_depth(params.max_depth),
                  try_split_models(model_spec::create_models(m_spec, 2)) {}

        unique_ptr<clustering_results>
        start(const unique_ptr<clustering_results> &cl_res,
              const clustering_input &input_params) override;

    private:
        cec_starter splitter;
        cross_entropy_clustering cec;
        const model_spec &m_spec;
        int max_k;
        int max_depth;
        vector<unique_ptr<model>> try_split_models;
        unique_ptr<clustering_results> try_split_cluster(const mat &x_mat);
    };
}

#endif //CEC_SPLIT_STARTER_H
