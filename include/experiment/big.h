#pragma once

#include "impl/message.h"
#include "impl/plaintext.h"

#include "util.h"

constexpr int LOGDELTA = 50;
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGQ = 1000; // max bit of a plaintext slot
constexpr int LOGN = 16;
constexpr int N = 1 << LOGN;
constexpr int K = 16;
constexpr int D = LOGN -1; // number of matrices multiplying

constexpr int L = 16; // number of RNS primes
constexpr uint64_t Q_primes[L] = { // 61bit primes, mod(q[i], 2^20) = 1
    0xffffffffffc0001ULL, 0x10000000006e0001ULL, 0xfffffffff840001ULL, 0x1000000000860001ULL,
    0xfffffffff6a0001ULL, 0x1000000000980001ULL, 0xfffffffff5a0001ULL, 0x1000000000b00001ULL,
    0x1000000000ce0001ULL, 0xfffffffff2a0001ULL, 0xfffffffff240001ULL, 0x1000000000f00001ULL,
    0xffffffffefe0001ULL, 0x10000000011a0001ULL, 0xffffffffeca0001ULL, 0xffffffffe9e0001ULL
};