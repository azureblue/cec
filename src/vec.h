#ifndef VEC_H
#define VEC_H

#include <algorithm>
#include <iostream>
#include <memory>

namespace cec {
    class mat;

    class row {
    public:
        const int size;

        row(row &&r) noexcept = default;

        row(double *ptr, int n)
                : size(n),
                  data_(ptr) {}

        inline double &operator[](int idx) {
            return data_[idx];
        }

        inline const double &operator[](int idx) const {
            return data_[idx];
        }

        void fill(double value) {
            std::fill(data_, data_ + size, value);
        }

        void set(const row &r) {
            this->operator=(r);
        }

        inline row &operator=(const row &r) {
            std::copy(r.data_, r.data_ + r.size, data_);
            return *this;
        }

        void operator*=(double value) {
            for (int i = 0; i < size; i++)
                data_[i] *= value;
        }

        void operator/=(double value) {
            for (int i = 0; i < size; i++)
                data_[i] /= value;
        }

        void operator+=(const row &r) {
            for (int i = 0; i < size; i++)
                data_[i] += r[i];
        }

        void operator-=(const row &r) {
            for (int i = 0; i < size; i++)
                data_[i] -= r[i];
        }

        double *data() {
            return data_;
        }

        const double *data() const {
            return data_;
        }

        friend std::ostream &operator<<(std::ostream &os, const row &r) {
            os << "{";
            for (int i = 0; i < r.size; i++)
                os << r[i] << (i == r.size - 1 ? "}" : ", ");
            return os;
        }

        static double dist(const row &a, const row &b) {
            double acc = 0;
            for (int i = 0; i < a.size; i++) {
                double diff = b[i] - a[i];
                acc += diff * diff;
            }
            return acc;
        }

    private:
        friend class mat;

        row(const row &r) = default;

        row sub(int offset, int n) {
            return {data_ + offset, n};
        }

        double *data_;
    };

    class mem {
    protected:
        explicit mem(int n)
                : ptr(new double[n]) {}

        mem(mem &&vm) = default;
        std::unique_ptr<double[]> ptr;
    };

    class vec: private mem, public row {
    public:
        explicit vec(int n)
                : mem(n),
                  row(ptr.get(), n) {}

        vec(const row &v)
                : vec(v.size) {
            (*this) = v;
        }

        vec(const vec &v)
                : vec(v.size) {
            (*this) = v;
        }

        vec(vec &&v) = default;

        using row::operator=;

        vec &operator=(const vec &v) {
            row::operator=(v);
            return *this;
        }
    };

    class mat {
    public:
        const int m, n;
    private:
        vec v;
    public:
        mat(int m, int n)
                : m(m),
                  n(n),
                  v(m * n) {}

        mat(const mat &ma)
                : mat(ma.m, ma.n) {
            (*this) = ma;
        }

        mat(mat &&ma) noexcept = default;

        inline const row operator[](int idx) const {
            return const_cast<mat *>(this)->operator[](idx);
        }

        inline row operator[](int idx) {
            return v.sub(n * idx, n);
        }

        inline mat &operator=(const mat &ma) {
            v.operator=(ma.v);
            return *this;
        }

        inline mat &operator=(mat &&ma) noexcept {
            *this = (ma);
            return *this;
        }

        double *data() {
            return v.data();
        }

        const double *data() const {
            return v.data();
        }

        void fill(double value) {
            v.fill(value);
        }

        void operator*=(double value) {
            v *= value;
        }

        void operator/=(double value) {
            v /= value;
        }

        void operator+=(const mat &m) {
            v += m.v;
        }

        void operator-=(const mat &m) {
            v -= m.v;
        }

        friend std::ostream &operator<<(std::ostream &os, const mat &m) {
            for (int j = 0; j < m.m; j++) {
                os << (j == 0 ? '(' : ' ');
                for (int k = 0; k < m.n; k++)
                    os << m[j][k] << (k < m.n - 1 ? ", " : "");
                os << (j == m.m - 1 ? ')' : '\n');
            }
            return os;
        }

        template<class M>
        class rows_iterator {
            friend class mat;

        public:
            typename std::conditional<std::is_const<M>::value, const row, row>::type
            operator*() {
                return ref[r];
            }

            void operator++() {
                r++;
            }

            bool operator==(const rows_iterator &ri) {
                return r == ri.r;
            }

            bool operator!=(const rows_iterator &ri) {
                return !operator==(ri);
            }

        private:
            rows_iterator(M &ref, int r)
                    : ref(ref),
                      r(r) {}

            M &ref;
            int r;
        };

        rows_iterator<mat> begin() {
            return {*this, 0};
        }

        rows_iterator<mat> end() {
            return {*this, m};
        }

        rows_iterator<const mat> begin() const {
            return {*this, 0};
        }

        rows_iterator<const mat> end() const {
            return {*this, m};
        }

        static mat outer_product(row &v) {
            int n = v.size;
            mat ma(n, n);
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++)
                    ma[j][k] = v[j] * v[k];
            return ma;
        }
    };
}

#endif /* VEC_H */

