#include <cstdint>

constexpr int LOGq = 60; // base modulus size
constexpr int LOGN = 16;
constexpr int N = 1 << LOGN;
constexpr int H = 192; // hamming weight of secret key

// for sanity check
constexpr int LOGDELTA = 50; // log delta used for mult
constexpr uint64_t Delta = 1ULL << LOGDELTA;
constexpr int LOGDELTA_boot_tilde = 34; // log delta for c2s
constexpr uint64_t Delta_boot_tilde = 1ULL << LOGDELTA_boot_tilde;