#ifndef CEC_RANDOM_H
#define CEC_RANDOM_H

#include <random>

namespace cec {
    class random {
    public:
        typedef std::mt19937 rand_gen;

        static rand_gen create_generator();

        static void set_seed(unsigned long seed);

    private:
        static unsigned long seed;
    };
}
#endif //CEC_RANDOM_H
