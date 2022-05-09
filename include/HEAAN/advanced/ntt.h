#pragma once

#include "HEAAN/Z_Q.h"
#include "mod.h"

#include <cstdint>

namespace {
inline int log(int N) {
    int i = 0, n = 1;
    while(n < N) {
        i++;
        n *= 2;
    }
    return i;
}

inline uint32_t bitReverse32(uint32_t x) {
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return ((x >> 16) | (x << 16));
}

inline uint32_t bitReverse(uint32_t x, int LOGN) {
    return bitReverse32(x) >> (32 - LOGN);
}
}

template <int N>
struct NTT {
    uint64_t q; // Object modulus, mod(q, 2N) = 1
    uint64_t psi_rev[N]; // psi be 2N-th root of q, psi^(reverse (0)) ~ psi^(reverse (N-1))
    uint64_t inv_psi_rev[N]; // psi_inv^(reverse (0)) ` psi_inv^(reverse(N-1))
    uint64_t inv_N; // N^-1 mod q

    NTT(uint64_t q, uint64_t psi);

    void ntt(uint64_t a[N]);
    void intt(uint64_t a[N]);
};

template<int N>
NTT<N>::NTT(uint64_t q, uint64_t psi) {
    int LOGN = log(N);
    uint64_t psi_orig[N], inv_psi_orig[N]; // power of psis in original order

    this->q = q;

    // compute psi^0 ~ psi^(N-1)
    psi_orig[0] = 1;
    for(int i = 1; i < N; ++i) {
        // psi_orig[i] = (psi * psi_orig[i-1]) % q
        psi_orig[i] = mul_mod(psi, psi_orig[i-1], q);
    }

    // set in reverse order
    for(int i = 0; i < N; ++i) {
        this->psi_rev[i] = psi_orig[bitReverse(i, LOGN)];
    }

    uint64_t inv_psi= inv_mod(psi, q);

    // compute inv_psi^0 ~ inv_psi^(N-1)
    inv_psi_orig[0] = 1;
    for(int i = 1; i < N; ++i) {
        inv_psi_orig[i] = mul_mod(inv_psi, inv_psi_orig[i-1], q);
    }

    // set in reverse order
    for(int i = 0; i < N; ++i) {
        this->inv_psi_rev[i] = inv_psi_orig[bitReverse(i, LOGN)];
    }

    this->inv_N = inv_mod(N, q);
}