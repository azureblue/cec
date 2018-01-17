#ifndef CEC_SPLIT_STARTER_H
#define CEC_SPLIT_STARTER_H

#include "starter.h"

namespace cec {
    struct split_starter_params {
        starter_params params;
        int starts;
        int split_tries;
        int max_k;
        int max_depth;
        int threads;
    };

    class split_starter {

    };
}

#endif //CEC_SPLIT_STARTER_H
