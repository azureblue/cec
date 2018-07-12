#ifndef CEC_COV_H
#define CEC_COV_H

#include "vec.h"

namespace cec {

    class mean: public vec {
    public:
        explicit mean(const mat &sample)
                : mean(sample.n) {
            for (auto &&p : sample) add_point(p);
            update();
        }

        explicit mean(int n)
                : vec(n),
                  acc(n) {
            acc.fill(0);
        }

        mean(const mean &initial) = default;

        mean(mean &&initial) noexcept = default;

        mean &operator=(const mean &m) = default;

        void add_point(const row &point) {
            acc += point;
            car++;
        }

        void rem_point(const row &point) {
            acc -= point;
            car--;
        }

        void update() {
            row::operator=(acc);
            (*this) /= car;
        }

        int card() const {
            return car;
        }

    private:
        using row::operator=;
        int car = 0;
        vec acc;
    };

    class covariance: public mat {
    public:
        covariance &operator=(const covariance &cov) = default;

        covariance(const covariance &) = default;

        covariance(covariance &&) noexcept = default;

        const cec::mean& mean() const {
            return mn;
        }

        int card() const {
            return mn.card();
        }

        static covariance estimate(const mat &sample) {
            cec::mean mn(sample);
            int n = sample.n;
            mat acc(n, n);
            acc.fill(0);
            vec t_vec(n);
            for (auto &&p : sample) {
                t_vec = p;
                t_vec -= mn;
                acc += mat::outer_product(t_vec);
            }
            acc /= sample.m;
            return covariance(acc, mn);
        }

        void add_point(const row &point) {
            int card = mn.card();
            double card_n = card + 1;
            cov_change(point, card / card_n, card / (card_n * card_n), *this);
            mn.add_point(point);
            mn.update();
        }

        void rem_point(const row &point) {
            int card = mn.card();
            double card_n = card - 1;
            cov_change(point, card / card_n, -card / (card_n * card_n), *this);
            mn.rem_point(point);
            mn.update();
        }

    protected:
        cec::mean mn;

        covariance(mat initial, cec::mean mn)
                : mat(std::move(initial)),
                  mn(std::move(mn)) {}

        inline void cov_change(const row &point, double cov_mul,
                               double new_cov_point_mul, mat &dst) {
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++)
                    dst[j][k] = (*this)[j][k] * cov_mul
                                + (mn[j] - point[j])
                                  * (mn[k] - point[k])
                                  * new_cov_point_mul;
        }
    };
}
#endif //CEC_COV_H
