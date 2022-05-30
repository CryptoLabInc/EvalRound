#pragma once

#include "impl/message.h"
#include "impl/simple_plaintext.h"
#include "impl/plaintext.h"

#include "setup.h"

#include <cstdint>

constexpr int LOGDELTA = 50;
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGQ = 1450; // max bit of a plaintext slot
constexpr int Diag = 16;