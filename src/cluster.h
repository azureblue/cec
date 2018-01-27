#ifndef CEC_CLUSTER_H
#define CEC_CLUSTER_H

#include "vec.h"
#include "models/model.h"

namespace cec {

    class mean : public vec {
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

    class covariance_mat : public mat {
    public:
        covariance_mat &operator=(const covariance_mat &cov) = default;

        covariance_mat(const covariance_mat &) = default;

        covariance_mat(covariance_mat &&) noexcept = default;

        const cec::mean &mean() const {
            return mn;
        }

        int card() const {
            return mn.card();
        }

        static covariance_mat estimate(const mat &sample) {
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
            return covariance_mat(acc, mn);
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

        covariance_mat(mat initial, cec::mean mn)
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

    class deferred_update_covariance : public covariance_mat {
    public:
        using mat::operator=;

        explicit deferred_update_covariance(const covariance_mat &initial)
                : covariance_mat(initial),
                  tmp_point(initial.m),
                  tmp_cov(initial) {}

        void add_point_tmp(const row &point) {
            tmp_point = point;
            card_change = 1;
            int card = covariance_mat::card();
            double card_n = card + 1;
            cov_change(point, card / card_n, card / (card_n * card_n), tmp_cov);
        }

        void rem_point_tmp(const row &point) {
            tmp_point = point;
            card_change = -1;
            int card = covariance_mat::card();
            double card_n = card - 1;
            cov_change(point, card / card_n, -card / (card_n * card_n), tmp_cov);
        }

        void apply_change() {
            (*this) = tmp_cov;
            if (card_change == 1) {
                mn.add_point(tmp_point);
                mn.update();
            } else if (card_change == -1){
                mn.rem_point(tmp_point);
                mn.update();
            }
            card_change = 0;
        }

        int tmp_card() const {
            return card() + card_change;
        }

        const mat &tmp_covariance() const {
            return tmp_cov;
        }

    private:
        int card_change = 0;
        vec tmp_point;
        mat tmp_cov;
    };

    class cluster {
    public:
        cluster(const model &mod, const covariance_mat &initial_covariance, const int m)
                : m(m),
                  mod(mod),
                  cov(initial_covariance),
                  eng(mod.energy(cov, cov.card(), m)),
                  tmp_eng(eng) {}

        double add_point(const row &point) {
            cov.add_point_tmp(point);
            tmp_eng = tmp_energy();
            return tmp_eng - eng;
        }

        double rem_point(const row &point) {
            cov.rem_point_tmp(point);
            tmp_eng = tmp_energy();
            return tmp_eng - eng;
        }

        const row &mean() const {
            return cov.mean();
        }

        int card() const {
            return cov.card();
        }

        const mat &covariance() const {
            return cov;
        }

        void apply_change() {
            cov.apply_change();
            eng = tmp_eng;
        }

        double energy() {
            return mod.energy(cov, cov.card(), m);
        }

    private:
        double tmp_energy() {
            return mod.energy(cov.tmp_covariance(), cov.tmp_card(), m);
        }

        const int m;
        const model &mod;
        deferred_update_covariance cov;
        double eng;
        double tmp_eng;
    };
}

#endif //CEC_CLUSTER_H
