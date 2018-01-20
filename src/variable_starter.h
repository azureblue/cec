#ifndef CEC_VARIABLE_STARTER_H
#define CEC_VARIABLE_STARTER_H

#include "cec_starter.h"

namespace cec {

    class variable_starter {
    public:

        variable_starter(vector<int> centers_number, unique_ptr<clustering_starter> cs)
                : centers_number(std::move(centers_number)),
                cs(std::move(cs)) {}

        unique_ptr<clustering_results>
        start(const mat &x, vector<shared_ptr<model_spec>> m_specs);

    private:
        vector<int> centers_number;
        unique_ptr<clustering_starter> cs;
    };
}

#endif //CEC_VARIABLE_STARTER_H
