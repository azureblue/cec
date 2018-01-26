#ifndef CEC_M_H
#define CEC_M_H

#include <cmath>

namespace cec {
    namespace m {
        const double PI = 3.14159265358979323846;
        const double E = 2.7182818284590452354;
        const double QNAN = std::numeric_limits<double>::quiet_NaN();
        const double INF = std::numeric_limits<double>::infinity();
        using std::log;
        using std::isnan;
    }
}
#endif //CEC_M_H
