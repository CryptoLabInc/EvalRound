#pragma once

#include <cstdint>
#include <cmath>

constexpr int LOGDELTA = 10;
constexpr std::uint64_t Delta = 1ULL << LOGDELTA;
constexpr int D = 1; // number of matrices multiplying
constexpr int LOGN = 10;
constexpr int N = 1 << LOGN;
constexpr int h = 64;
constexpr int K = 16;
