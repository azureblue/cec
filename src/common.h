#ifndef CEC_STD_COMMON_H
#define CEC_STD_COMMON_H

#include <vector>
#include <string>
#include <memory>

namespace cec {
    using std::vector;
    using std::string;
    using std::unique_ptr;
    using std::shared_ptr;
    using std::make_shared;

#if defined(__cplusplus) && __cplusplus >= 201402L || defined(__cpp_lib_make_unique)
    using std::make_unique
#else
    template <typename T, typename ...Args>
    std::unique_ptr<T> make_unique(Args&& ...args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
#endif
}

#endif //CEC_STD_COMMON_H
