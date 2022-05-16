#pragma once

#include "impl/message.h"
#include "impl/plaintext.h"

#include "setup.h"

constexpr int LOGDELTA = 60; // log delta used for operations except c2s
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGDELTA_TILDE = 30; // log delta only for c2s
constexpr uint64_t Delta_tilde = 1ULL << LOGDELTA;
constexpr int LOGQ = 1000; // max bit of a plaintext slot
constexpr int LOGN = 9;
constexpr int N = 1 << LOGN;

constexpr int K = 24; // max number of periods of sine function to approximate on EvalMod
constexpr int H = 128; // hamming weight of secret key