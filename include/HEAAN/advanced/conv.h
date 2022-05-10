#pragma once

#include "HEAAN/advanced/crt.h"
#include "HEAAN/advanced/ntt.h"

#include <memory>

template<int L, int N>
struct CONV{
    CRT<L> crt;
    std::unique_ptr<NTT<N>> ntt[L];

    CONV(const uint64_t q[L]);
};

template<int L, int N>
CONV<L, N>::CONV(const uint64_t q[L]) : crt{q} {
    for(int i = 0; i < L; ++i) {
        ntt[i] = std::make_unique<NTT<N>>(q[i]);
    }
}

