#pragma once

#include "impl/message.h"
#include "impl/plaintext.h"

#include "setup.h"

#include <cstdint>
    
constexpr int LOGq = 60; // base modulus size
constexpr int LOGQ = 1450; // max bit of a plaintext slot
constexpr int LOGN = 17;
constexpr int N = 1 << LOGN;

constexpr int LOGDELTA = 50; // log delta used for mult
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGDELTA_boot = LOGq; // log delta for bootstrap
constexpr uint64_t Delta_boot = 1ULL << LOGDELTA_boot;
constexpr int LOGDELTA_boot_tilde = 28; // log delta for c2s
constexpr uint64_t Delta_boot_tilde = 1ULL << LOGDELTA_boot_tilde;

constexpr int K = 24; // max number of periods of sine function to approximate on EvalMod
constexpr int H = 128; // hamming weight of secret key
constexpr int G = 4; // number of grouped multiplication of encoding matrix
