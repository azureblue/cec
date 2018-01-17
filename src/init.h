#ifndef CEC_INIT_H
#define CEC_INIT_H

#include <utility>

#include "vec.h"
#include "random.h"
#include "common.h"

namespace cec {

    using std::vector;
    using std::unique_ptr;
    using std::shared_ptr;

    class centers_init {
    public:
        virtual ~centers_init() = default;

        virtual mat init(const mat &x, int k) = 0;
    };

    class random_init : public centers_init {
    public:
        random_init()
                : gen(random::create_generator()) {}

        mat init(const mat &x, int k) override;

    private:
        random::rand_gen gen;
    };

    class fixed_init : public centers_init {
    public:
        explicit fixed_init(mat c_mat)
                : c_mat(std::move(c_mat)) {}

        mat init(const mat &x, int k) override;

    private:
        mat c_mat;
    };

    class kmeanspp_init : public centers_init {
    public:
        mat init(const mat &x, int k) override;

        explicit kmeanspp_init()
                : gen(random::create_generator()) {}

    private:
        random::rand_gen gen;
        vector<double> dists;
        vector<double> sums;
    };

    class centers_init_spec {
    public:
        virtual ~centers_init_spec() = default;

        virtual unique_ptr<centers_init> create() const = 0;
    };

    class random_init_spec : public centers_init_spec {
    public:
        unique_ptr<centers_init> create() const override;
    };

    class fixed_init_spec : public centers_init_spec {
    public:
        const mat c_mat;

        explicit fixed_init_spec(mat c_mat)
                : c_mat(std::move(c_mat)) {}

        unique_ptr<centers_init> create() const override;
    };

    class kmeanspp_init_spec : public centers_init_spec {
    public:
        unique_ptr<centers_init> create() const override;
    };

    class assignment_init {
    public:
        virtual ~assignment_init() = default;

        virtual std::vector<int> init(const mat &x, const mat &c) = 0;
    };

    class closest_assignment : public assignment_init {
    public:
        std::vector<int> init(const mat &x, const mat &c) override;
    };
}

#endif //CEC_INIT_H
