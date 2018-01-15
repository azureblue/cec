#ifndef CEC_RANDOM_H
#define CEC_RANDOM_H

#include <random>

namespace cec {
    class random {
    public:
        typedef std::mt19937 rand_gen;
        typedef rand_gen::result_type result_type;

        static rand_gen create_generator();

        static void set_seed(result_type seed) noexcept;

    private:
        static unsigned long seed;
    };
}
#endif //CEC_RANDOM_H
