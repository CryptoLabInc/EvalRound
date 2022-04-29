#pragma once

#include "util.h"

#include <cstdint>
#include <cmath>

constexpr int LOGDELTA = 10;
constexpr std::uint64_t Delta = 1ULL << LOGDELTA;
constexpr std::uint64_t DeltaSq = 1ULL << (2 * LOGDELTA);
constexpr std::uint64_t DeltaTr = 1ULL << (3 * LOGDELTA);
constexpr int LOGQ = 24;
constexpr int LOGQ_UP = LOGQ * 2;
constexpr int N = 1 << 10;
constexpr int h = 64;
constexpr int K = 16;