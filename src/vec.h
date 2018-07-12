#ifndef VEC_H
#define VEC_H

#include <algorithm>
#include <iostream>
#include <memory>

namespace cec {
    class vec;
    class mat;
    class row {
    public:
        const int size;

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

        static double dist_sq(const row &a, const row &b) {
            double acc = 0;
            for (int i = 0; i < a.size; i++) {
                double diff = b[i] - a[i];
                acc += diff * diff;
            }
            return acc;
        }

    private:
        friend class vec;
        friend class mat;

        row(const row &r) = default;

        row(row &&r) noexcept = default;

        row(double *ptr, int n)
                : size(n),
                  data_(ptr) {}


        row sub(int offset, int n) {
            return {data_ + offset, n};
        }

        double *data_;
    };

    template<typename T>
    class storage {
    protected:
        explicit storage(int n)
                : ptr(new T[n]) {}

        storage(storage &&vm) noexcept = default;

        inline T *get_storage() {
            return ptr.get();
        }
    private:
        std::unique_ptr<T[]> ptr;
    };

    class vec: private storage<double>, public row {
    public:
        explicit vec(int n)
                : storage<double>(n),
                  row(get_storage(), n) {}

        explicit vec(const row &v)
                : vec(v.size) {
            operator=(v);
        }

        vec(const vec &v)
                : vec(v.size) {
            operator=(v);
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
        vec data_vec;
    public:
        mat(int m, int n)
                : m(m),
                  n(n),
                  data_vec(m * n) {}

        mat(const mat &ma)
                : mat(ma.m, ma.n) {
            (*this) = ma;
        }

        mat(mat &&ma) noexcept = default;

        inline const row operator[](int idx) const {
            return const_cast<mat *>(this)->operator[](idx);
        }

        inline row operator[](int idx) {
            return data_vec.sub(n * idx, n);
        }

        inline mat &operator=(const mat &ma) {
            data_vec.operator=(ma.data_vec);
            return *this;
        }

        inline mat &operator=(mat &&ma) noexcept {
            *this = ma;
            return *this;
        }

        double *data() {
            return data_vec.data();
        }

        const double *data() const {
            return data_vec.data();
        }

        void fill(double value) {
            data_vec.fill(value);
        }

        void operator*=(double value) {
            data_vec *= value;
        }

        void operator/=(double value) {
            data_vec /= value;
        }

        void operator+=(const mat &m) {
            data_vec += m.data_vec;
        }

        void operator-=(const mat &m) {
            data_vec -= m.data_vec;
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
            using difference_type = long;
            using value_type = typename std::conditional<std::is_const<M>::value, const row, row>::type;
            using pointer = value_type*;
            using reference = value_type&;
            using iterator_category = std::forward_iterator_tag;

            value_type
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
