#pragma once

constexpr int LOGDELTA = 50;
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGQ = 512;
constexpr int LOGN = 16;
constexpr int N = 1 << LOGN;
constexpr int K = 16;
constexpr int D = LOGN -1; // number of matrices multiplying