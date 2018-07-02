#ifndef R_EXT_PTR_H
#define R_EXT_PTR_H

#include <Rinternals.h>
#include <Rdefines.h>

namespace cec {
    namespace r {
        template<typename T>
        class r_ext_ptr {
        public:
            r_ext_ptr() {
                r_ptr = PROTECT(R_MakeExternalPtr(nullptr, NULL, NULL));
                R_RegisterCFinalizerEx(r_ptr, r_ext_ptr::finalize, TRUE);
            }

            r_ext_ptr(r_ext_ptr<T> &r_p) = delete;

            r_ext_ptr(r_ext_ptr<T> &&r_p) noexcept {
                (*this) = std::move(r_p);
            }

            bool operator==(const r_ext_ptr &r_p) const {
                return t_ptr == r_p.t_ptr;
            }

            bool operator!=(const r_ext_ptr &r_p) const {
                return !(r_p == *this);
            }

            explicit operator bool() {
                return t_ptr != nullptr;
            }

            r_ext_ptr &operator=(r_ext_ptr &r_p) = delete;

            r_ext_ptr &operator=(r_ext_ptr &&r_p) noexcept {
                r_ptr = r_p.r_ptr;
                t_ptr = r_p.t_ptr;
                r_p.r_ptr = NULL;
                r_p.t_ptr = nullptr;
                return *this;
            }

            virtual ~r_ext_ptr() {
                if (r_ptr == NULL)
                    return;
                finalize(r_ptr);
                UNPROTECT_PTR(r_ptr);
            }

            template <typename ...Args>
            void init(Args &&...args) {
                finalize(r_ptr);
                t_ptr = new T(std::forward<Args>(args)...);
                R_SetExternalPtrAddr(r_ptr, t_ptr);
            }

            void reset(T *ptr) {
                finalize(r_ptr);
                t_ptr = ptr;
                R_SetExternalPtrAddr(r_ptr, ptr);
            }

            const T *get() const {
                return t_ptr;
            }

            const T *operator->() const {
                return t_ptr;
            }

            const T &operator*() const {
                return *t_ptr;
            }

            T *get() {
                return t_ptr;
            }

            T *operator->() {
                return t_ptr;
            }

            T &operator*() {
                return *t_ptr;
            }

        private:
            T *t_ptr = nullptr;
            SEXP r_ptr = NULL;

            static void finalize(SEXP r_ptr) {
                T *ptr = (T *) R_ExternalPtrAddr(r_ptr);
                if (ptr == nullptr)
                    return;
                delete ptr;
                R_ClearExternalPtr(r_ptr);
            }
        };

        template<typename T, typename ...Args>
        static r_ext_ptr<T> make_r_ext(Args &&... args) {
            r_ext_ptr<T> ptr;
            ptr.init(std::forward<Args>(args)...);
            return ptr;
        }

    }
}

#endif //R_EXT_PTR_H
