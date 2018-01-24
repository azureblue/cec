#ifndef CEC_VARIABLE_STARTER_H
#define CEC_VARIABLE_STARTER_H

#include "cec_starter.h"
#include "parallel_starter.h"

namespace cec {
    class variable_starter {
    public:
        using clustering_function = std::function<unique_ptr<clustering_results>
                (const mat&, const vector<shared_ptr<model_spec>>&)>;

        variable_starter(clustering_function &&cl_starter, vector<int> centers_number)
                : cl_starter(std::move(cl_starter)),
                  centers_number(centers_number) {}

        unique_ptr<clustering_results>
        start(const mat &x, vector<shared_ptr<model_spec>> m_specs);

    private:
        clustering_function cl_starter;
        vector<int> centers_number;
    };
}

#endif //CEC_VARIABLE_STARTER_H
