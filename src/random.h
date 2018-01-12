#ifndef CEC_RANDOM_H
#define CEC_RANDOM_H

#include <random>

namespace cec {
    namespace random {
        class context {
        public:
            static std::mt19937 create_generator();
            static void set_seed(unsigned long seed);

        private:
            static unsigned long seed;
        };
    }
}
#endif //CEC_RANDOM_H
