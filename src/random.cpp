#include "random.h"
unsigned long cec::random::context::seed = 0;

std::mt19937 cec::random::context::create_generator() {
    std::mt19937 gen(context::seed);
    context::seed++;
    return gen;
}

void cec::random::context::set_seed(unsigned long seed) {
    context::seed = seed;
}
