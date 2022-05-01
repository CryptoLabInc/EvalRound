#include <random>

namespace {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1.0, 1.0);
}

template <class Iter>
void sampleUniform(Iter begin, Iter end) {
    for(Iter it = begin; it != end; ++it) {
        *it = dist(gen);
    }
}