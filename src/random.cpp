#include "random.h"
unsigned long cec::random::seed = 0;

cec::random::rand_gen cec::random::create_generator() {
    std::mt19937 mt(seed);
    seed++;
    return mt;
}

void cec::random::set_seed(result_type seed) noexcept {
    random::seed = seed;
}
