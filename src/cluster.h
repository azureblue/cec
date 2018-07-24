#ifndef CEC_CLUSTER_H
#define CEC_CLUSTER_H

#include "vec.h"
#include "models/model.h"

namespace cec {
    class deferred_update_covariance : public covariance {
    public:
        using covariance::operator=;

        explicit deferred_update_covariance(const covariance &initial)
                : covariance(initial),
                  tmp_point(initial.m),
                  tmp_cov(initial) {}

        void add_point_tmp(const row &point) {
            tmp_cov = *this;
            tmp_cov.add_point(point);
        }

        void rem_point_tmp(const row &point) {
            tmp_cov = *this;
            tmp_cov.rem_point(point);
        }

        void apply_change() {
            (*this) = tmp_cov;
        }

        const covariance &tmp_covariance() const {
            return tmp_cov;
        }

    private:
        vec tmp_point;
        covariance tmp_cov;
    };

    class cluster {
    public:
        cluster(const model &mod, const covariance &initial_covariance, const int m)
                : m(m),
                  mod(mod),
                  cov(initial_covariance),
                  eng(mod.energy(cov, m)),
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
            return mod.energy(cov, m);
        }

    private:
        double tmp_energy() {
            return mod.energy(cov.tmp_covariance(), m);
        }

        const int m;
        const model &mod;
        deferred_update_covariance cov;
        double eng;
        double tmp_eng;
    };
}

#endif //CEC_CLUSTER_H
