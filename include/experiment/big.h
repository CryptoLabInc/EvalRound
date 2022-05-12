#pragma once

#include "impl/message.h"
#include "impl/plaintext.h"

#include "util.h"

constexpr int LOGDELTA = 50;
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGQ = 250; // max bit of a plaintext slot
constexpr int LOGN = 16;
constexpr int N = 1 << LOGN;
constexpr int K = 16;
constexpr int D = LOGN -1; // number of matrices multiplying