#pragma once

#include <cstdint>
#include <cmath>

constexpr int LOGDELTA = 10;
constexpr std::uint64_t Delta = 1ULL << LOGDELTA;
constexpr std::uint64_t DeltaSq = 1ULL << (2 * LOGDELTA);
constexpr int D = 1; // number of matrices multiplying
constexpr std::uint64_t DeltaTotal = 1ULL << (LOGDELTA * (D + 1));
constexpr int LOGQ = 500;
constexpr int LOGQ_UP = LOGQ * 2;
constexpr int LOGN = 10;
constexpr int N = 1 << LOGN;
constexpr int h = 64;
constexpr int K = 16;
