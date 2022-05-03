#pragma once

#include <cstdint>
#include <cmath>

constexpr int LOGDELTA = 20;
constexpr std::uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGN = 16;
constexpr int N = 1 << LOGN;
constexpr int K = 16;
