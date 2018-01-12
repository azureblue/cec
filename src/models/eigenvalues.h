#ifndef EIGENVALUES_H
#define EIGENVALUES_H

#include "cov_utils.h"
#include "model.h"

namespace cec {
    class eigenvalues : public model {
    public:
        explicit eigenvalues(int n, std::vector<double> values)
                : n(n),
                  given_values(std::move(values)),
                  eigenvalues_calc(n),
                  tmp_values(n),
                  ce_constant(std::log(std::pow(2.0 * m::PI, n)
                                    * product(eigenvalues::given_values)) / 2.0) {}

        double cross_entropy(const mat &cov) const noexcept override {
            if (!eigenvalues_calc.eigenvalues(cov, tmp_values.data()))
                return m::QNAN;
            double values_ratio_sum = 0;
            for (int i = 0; i < n; i++)
                values_ratio_sum += tmp_values[i] / given_values[i];
            return ce_constant + values_ratio_sum / 2.0;
        }

    private:
        const int n;
        const std::vector<double> given_values;
        const eigenvalues_calculator eigenvalues_calc;
        mutable std::vector<double> tmp_values;
        const double ce_constant;

        static double product(const std::vector<double> &values) {
            double prod = 1;
            for (auto &&v : values)
                prod *= v;
            return prod;
        }
    };
}
#endif /* EIGENVALUES_H */

