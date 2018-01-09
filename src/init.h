#ifndef CEC_INIT_H
#define CEC_INIT_H

#include "vec.h"

namespace cec {
    class centers_init {
    public:
        virtual mat init(const mat &x, int k) = 0;
    };

    class random_init : public centers_init {
    public:
        mat init(const mat &x, int k) override;

        random_init()
                : mt(std::random_device()()) {}

    private:
        std::mt19937 mt;
    };

    class kmeanspp_init : public centers_init {
    public:
        mat init(const mat &x, int k) override;

        kmeanspp_init(int m_max)
                : mt(std::random_device()()),
                  dists(m_max),
                  sums(m_max) {}

    private:
        std::mt19937 mt;
        std::vector<double> dists;
        std::vector<double> sums;
    };

    class assignment_init {
    public:
        virtual std::vector<int> init(const mat &x, const mat &c) = 0;
    };

    class closest_assignment : public assignment_init {
    public:
        std::vector<int> init(const mat &x, const mat &c) override;
    };
}

#endif //CEC_INIT_H
