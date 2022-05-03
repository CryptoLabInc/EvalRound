#pragma once

#include <cstdint>
#include <cmath>

constexpr int LOGDELTA = 50;
constexpr std::uint64_t Delta = 1ULL << LOGDELTA;
constexpr int D = 3; // number of matrices multiplying
constexpr int LOGQ = 256;
constexpr int LOGQ_UP = LOGQ * 2;
constexpr int LOGN = 10;
constexpr int N = 1 << LOGN;
constexpr int h = 64;
constexpr int K = 16;
