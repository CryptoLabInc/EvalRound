#pragma once

#include <cstdint>
#include <cmath>

constexpr std::uint64_t Delta = 1ULL << 50;
constexpr int LOGQ = 500;
constexpr int N = 1 << 3;

void set_test_message(double zr[], double zi[]) {
    for(int i = 0; i < N / 2; ++i) {
        double x = (double) i / (N/2) * M_PI;
        zr[i] = cos(x);
        zi[i] = sin(x);
    }
}