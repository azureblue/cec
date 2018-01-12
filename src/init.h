#ifndef CEC_INIT_H
#define CEC_INIT_H

#include <utility>

#include "vec.h"
#include "random.h"

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
                : mt(random::context::create_generator()) {}

    private:
        std::mt19937 mt;
    };

    class fixed_init : public centers_init {
    public:
        mat init(const mat &x, int k) override;

        fixed_init(mat c_mat)
                : c_mat(c_mat) {}

    private:
        mat c_mat;
    };

    class kmeanspp_init : public centers_init {
    public:
        mat init(const mat &x, int k) override;

        kmeanspp_init(int m_max)
                : mt(random::context::create_generator()),
                  dists(m_max),
                  sums(m_max) {}

    private:
        std::mt19937 mt;
        std::vector<double> dists;
        std::vector<double> sums;
    };

    class centers_init_spec {
    public:
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
    public:
        virtual unique_ptr<assignment_init> create() = 0;
    };

    class closest_init_spec: public assignment_init_spec {
        unique_ptr<assignment_init> create() override {
            return unique_ptr<assignment_init>(new closest_assignment);
        }
    };

    class random_init_spec: public centers_init_spec {
    public:
        unique_ptr<centers_init> create() override {
            return unique_ptr<centers_init>(new random_init);
        }
    };

    class kmeanspp_init_spec: public centers_init_spec {
    public:
        const int max_m;

        kmeanspp_init_spec(const int max_m)
                : max_m(max_m) {}

        unique_ptr<centers_init> create() override {
            return unique_ptr<centers_init>(new kmeanspp_init(max_m));
        }
    };

    class initializer {
    public:
        initializer(unique_ptr<centers_init> &&c_init, unique_ptr<assignment_init> &&a_init)
                : c_init(std::move(c_init)),
                  a_init(std::move(a_init)) {}

        vector<int> init(const mat &x, int k) {
            return a_init->init(x, c_init->init(x, k));
        }

    private:
        unique_ptr<centers_init> c_init;
        unique_ptr<assignment_init> a_init;
    };

    class init_spec {
    public:
        init_spec(shared_ptr<centers_init_spec> ci, shared_ptr<assignment_init_spec> ai)
                : ci(std::move(ci)),
                  ai(std::move(ai)) {}

        unique_ptr<initializer> create() {
            return unique_ptr<initializer>(new initializer(ci->create(), ai->create()));
        }

        shared_ptr<centers_init_spec> ci;
        shared_ptr<assignment_init_spec> ai;
    };
}

#endif //CEC_INIT_H
