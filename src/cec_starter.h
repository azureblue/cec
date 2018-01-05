#ifndef CEC_STARTER_H
#define CEC_STARTER_H

#include "vec.h"
#include "single_start_input.h"

namespace cec {

    class cec_starter {
    public:
        cec_starter() = default;

        single_start_results start(const single_start_input &in);
        std::vector<mat> split_points(const mat &points, const std::vector<int> &assignment, int k);
    };


}

#endif //CEC_STARTER_H