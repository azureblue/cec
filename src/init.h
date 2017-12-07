#ifndef CEC_INIT_H
#define CEC_INIT_H

#include "vec.h"

namespace cec {
    class centers_init {
    public:
        virtual mat init(mat x, int k) = 0;
    };

    class random_init : public centers_init {
    public:
        mat init(mat x_mat, int k) override {
            mat c_mat(k, x_mat.n);
            std::uniform_int_distribution<int> unif_int(0, x_mat.m - 1);
            for (int i = 0; i < k; i++)
                c_mat[i] = x_mat[unif_int(mt)];
            return c_mat;
        }

        random_init()
                : rd(),
                  mt(rd()) {}
    private:
        std::random_device rd;
        std::mt19937 mt;

    };

    class assignment_init {
    public:
        virtual std::vector<int> init(const mat &x_mat, const mat &c_mat) = 0;
    };

    class closest_assignment : public assignment_init {
    public:
        std::vector<int> init(const mat &x_mat, const mat &c_mat) override {
            int m = x_mat.m;
            int k = c_mat.m;
            std::vector<int> asgn(m);
            for (int i = 0; i < m; i++) {
                double b_dist = std::numeric_limits<double>::max();
                const vec &point = x_mat[i];
                int b_row = -1;
                for (int j = 0; j < k; j++) {
                    double dist = vec::dist(point, c_mat[j]);
                    if (dist < b_dist) {
                        b_dist = dist;
                        b_row = j;
                    }
                }
                asgn[i] = b_row;
            }
            return asgn;
        }
    };
};

#endif //CEC_INIT_H
