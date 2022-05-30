#pragma once

#include "HEAAN/arith/Z_Q.h"
#include "util/util.h"
#include "mod.h"

namespace {
void findPrimeFactors(std::set<uint64_t> &s, uint64_t n) {
    s.clear();

    while (n % 2 == 0) {
        s.insert(2);
        n /= 2;
    }

    for (uint64_t i = 3; i * i <= n; i += 2) {
        while (n % i == 0) {
            s.insert(i);
            n /= i;
        }
    }

    if (n > 2)
        s.insert(n);
}

// find primitive root for q
uint64_t findPrimitiveRoot(uint64_t prime) {
    std::set<uint64_t> s;
    uint64_t phi = prime - 1;
    findPrimeFactors(s, phi);
    for (uint64_t r = 2; r <= phi; r++) {
        bool passed = true;
        for (uint64_t it : s) {
            if (pow_mod(r, phi / it, prime) == 1) {
                passed = false;
                break;
            }
        }

        if (passed)
            return r;
    }

    return 0; // failed to find
}
}

template <int N>
struct NTT {
    uint64_t q; // Object modulus, mod(q, 2N) = 1
    uint64_t psi_rev[N]; // psi be 2N-th root of q, psi^(reverse (0)) ~ psi^(reverse (N-1))
    uint64_t inv_psi_rev[N]; // psi_inv^(reverse (0)) ` psi_inv^(reverse(N-1))
    uint64_t inv_N; // N^-1 mod q

    NTT(uint64_t q);

    void ntt(uint64_t a[N]) const;
    void intt(uint64_t a[N]) const;
};

template<int N>
NTT<N>::NTT(uint64_t q) {
    int LOGN = log(N);
    uint64_t psi_orig[N], inv_psi_orig[N]; // power of psis in original order

    this->q = q;

    // compute psi
    uint64_t psi = findPrimitiveRoot(q);
    psi = pow_mod(psi, (q-1) / (2*N), q);

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

template<int N>
void NTT<N>::ntt(uint64_t a[N]) const {
    for(int m = 1, t = N/2; m < N; m*=2, t/=2) {
        for(int i =0; i <m; i++) {
            int j1 = 2*i*t;
            int j2 = j1 + t - 1;
            uint64_t S = psi_rev[m+i];
            for(int j = j1; j <= j2; j++) {
                uint64_t U = a[j];
                uint64_t V = mul_mod(a[j+t], S, q);
                a[j] = add_mod(U, V, q);
                a[j+t] = sub_mod(U, V, q);
            }
        }
    }
}

template<int N>
void NTT<N>::intt(uint64_t a[N]) const {
    for(int t = 1, m = N; m > 1; m /=2, t *= 2) {
        int j1 = 0;
        int h = m/2;
        for(int i = 0; i < h; ++i) {
            int j2 = j1 + t - 1;
            uint64_t S = inv_psi_rev[h+i];
            for(int j = j1; j <= j2; ++j) {
                uint64_t U = a[j];
                uint64_t V = a[j+t];
                a[j] = add_mod(U, V, q);
                a[j+t] = mul_mod(sub_mod(U, V, q), S, q);
            }
            j1 += 2*t;
        }
    }
    for(int j =0; j < N; ++j) {
        a[j] = mul_mod(a[j], inv_N, q);
    }
}