#pragma once

#include "impl/message.h"
#include "impl/plaintext.h"

#include "setup.h"

#include <cstdint>

constexpr int LOGDELTA = 50;
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGQ = 1000; // max bit of a plaintext slot
constexpr int LOGN = 15;
constexpr int N = 1 << LOGN;
