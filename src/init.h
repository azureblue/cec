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
        virtual ~centers_init() = default;
        virtual mat init(const mat &x, int k) = 0;
    };

    class random_init : public centers_init {
    public:
        random_init();
        mat init(const mat &x, int k) override;

    private:
        std::mt19937 mt;
    };

    class fixed_init : public centers_init {
    public:
        explicit fixed_init(mat c_mat);
        mat init(const mat &x, int k) override;

    private:
        mat c_mat;
    };

    class kmeanspp_init : public centers_init {
    public:
        mat init(const mat &x, int k) override;
        explicit kmeanspp_init();

    private:
        std::mt19937 mt;
        std::vector<double> dists;
        std::vector<double> sums;
    };

    class centers_init_spec {
    public:
        virtual ~centers_init_spec() = default;
        virtual unique_ptr<centers_init> create() = 0;
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

    class assignment_init_spec {
    public:
        virtual ~assignment_init_spec() = default;
        virtual unique_ptr<assignment_init> create() = 0;
    };

    class closest_init_spec: public assignment_init_spec {
        unique_ptr<assignment_init> create() override;
    };

    class random_init_spec: public centers_init_spec {
    public:
        unique_ptr<centers_init> create() override;
    };

    class fixed_init_spec: public centers_init_spec {
    public:
        const mat c_mat;

        explicit fixed_init_spec(const mat &c_mat);

        unique_ptr<centers_init> create() override;
    };

    class kmeanspp_init_spec: public centers_init_spec {
    public:
        unique_ptr<centers_init> create() override;
    };

    class initializer {
    public:
        initializer(unique_ptr<centers_init> &&c_init, unique_ptr<assignment_init> &&a_init);

        vector<int> init(const mat &x, int k);

    private:
        unique_ptr<centers_init> c_init;
        unique_ptr<assignment_init> a_init;
    };

    class initializer_spec {
    public:
        initializer_spec(shared_ptr<centers_init_spec> ci, shared_ptr<assignment_init_spec> ai);

        unique_ptr<initializer> create() const;

        shared_ptr<centers_init_spec> ci;
        shared_ptr<assignment_init_spec> ai;
    };
}

#endif //CEC_INIT_H
