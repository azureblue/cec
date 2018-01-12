#ifndef CEC_INIT_H
#define CEC_INIT_H

#include "vec.h"

namespace cec {

    using std::vector;
    using std::unique_ptr;
    using std::shared_ptr;

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

    class centers_init_spec {
        virtual unique_ptr<centers_init> create() = 0;
    };



    class assignment_init {
    public:
        virtual std::vector<int> init(const mat &x, const mat &c) = 0;
    };

    class closest_assignment : public assignment_init {
    public:
        std::vector<int> init(const mat &x, const mat &c) override;
    };

    class assignment_init_spec {
        virtual unique_ptr<assignment_init> create() = 0;
    };

    class closest_init_spec: public assignment_init_spec {
        unique_ptr<assignment_init> create() override {
            return unique_ptr<assignment_init>(new closest_assignment);
        }
    };
}

#endif //CEC_INIT_H
