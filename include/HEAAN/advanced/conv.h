#pragma once

#include "crt.h"
#include "ntt.h"
#include "mod.h"

#include <memory>

template<int L, int LOGQ, int N>
struct CONV{
    CRT<L, LOGQ> crt;
    std::unique_ptr<NTT<N>> ntt[L];

    CONV(const uint64_t q[L]);

    void conv( const R_Q<LOGQ,N>& A,
	       const R_Q<LOGQ,N>& B, R_Q<LOGQ,N>& C ) const;

    private:
        void icrt_polynomial(const R_Q<LOGQ, N>& A, uint64_t A_rns[L][N]) const;
        void crt_polynomial(const uint64_t A_rns[L][N], R_Q<LOGQ, N> &A) const;
};

template<int L, int LOGQ, int N>
CONV<L, LOGQ, N>::CONV(const uint64_t q[L]) : crt{q} {
    for(int i = 0; i < L; ++i) {
        ntt[i] = std::make_unique<NTT<N>>(q[i]);
    }
}

template<int L, int LOGQ, int N>
void CONV<L, LOGQ, N>::conv( const R_Q<LOGQ,N>& A,
	       const R_Q<LOGQ,N>& B, R_Q<LOGQ,N>& C ) const {
    uint64_t A_rns[L][N], B_rns[L][N], C_rns[L][N];

    icrt_polynomial(A, A_rns);
    icrt_polynomial(B, B_rns);

    for(int i = 0; i < L; ++i) {
        ntt[i]->ntt(A_rns[i]);
        ntt[i]->ntt(B_rns[i]);
        for(int j = 0; j < N; ++j) {
            C_rns[i][j] = mul_mod(A_rns[i][j], B_rns[i][j] , ntt[i]->q);
        }
        ntt[i]->intt(C_rns[i]);
    }

    crt_polynomial(C_rns, C);
}

template<int L, int LOGQ, int N>
void CONV<L, LOGQ, N>::icrt_polynomial(const R_Q<LOGQ, N>& A, uint64_t A_rns[L][N]) const {
    for(int i = 0; i < N; ++i) {
        uint64_t Ai_rns[L];
        crt.icrt(A[i], Ai_rns);
        for(int j = 0; j < L; ++j){
            A_rns[j][i] = Ai_rns[j];
        }
    }
}

template<int L, int LOGQ, int N>
void CONV<L, LOGQ, N>::crt_polynomial(const uint64_t A_rns[L][N], R_Q<LOGQ, N> &A) const {
    for(int i = 0; i < N; ++i) {
        uint64_t Ai_rns[L];
        for(int j = 0; j < L; ++j) {
            Ai_rns[j] = A_rns[j][i];
        }
        crt.crt(Ai_rns, A[i]);
    }
}



template< int LOGQ, int N >
void conv( const R_Q<LOGQ,N>& A,
	       const R_Q<LOGQ,N>& B, R_Q<LOGQ,N>& C ){
    // we fix params for primes since this will be used only for conv.
    // 102 60bit primes, mod(q[i], 2^18) = 1
    static constexpr uint64_t Q_primes[] = {
        0xffffffffffc0001ULL,
        0xfffffffff840001ULL,
        0x1000000000980001ULL,
        0x1000000000b00001ULL,
        0xfffffffff240001ULL,
        0x1000000000f00001ULL,
        0xffffffffe7c0001ULL,
        0xffffffffe740001ULL,
        0x1000000001a00001ULL,
        0xffffffffe4c0001ULL,
        0xffffffffe440001ULL,
        0xffffffffe400001ULL,
        0x1000000002340001ULL,
        0xffffffffdbc0001ULL,
        0xffffffffd840001ULL,
        0x1000000002940001ULL,
        0xffffffffd680001ULL,
        0xffffffffd000001ULL,
        0xffffffffcf00001ULL,
        0xffffffffcdc0001ULL,
        0xffffffffcc40001ULL,
        0x1000000003680001ULL,
        0x1000000003900001ULL,
        0xffffffffc300001ULL,
        0x1000000003ec0001ULL,
        0xffffffffbf40001ULL,
        0xffffffffbdc0001ULL,
        0xffffffffb880001ULL,
        0x1000000004d40001ULL,
        0x1000000004f80001ULL,
        0x1000000005040001ULL,
        0xffffffffaec0001ULL,
        0x1000000005400001ULL,
        0x1000000005480001ULL,
        0x10000000054c0001ULL,
        0x1000000005640001ULL,
        0xffffffffa380001ULL,
        0xffffffffa200001ULL,
        0xffffffffa0c0001ULL,
        0x10000000069c0001ULL,
        0xffffffff9600001ULL,
        0x1000000006e00001ULL,
        0xffffffff91c0001ULL,
        0xffffffff8f40001ULL,
        0x10000000075c0001ULL,
        0x1000000007640001ULL,
        0xffffffff8680001ULL,
        0x1000000007d40001ULL,
        0x10000000080c0001ULL,
        0xffffffff7e40001ULL,
        0xffffffff7bc0001ULL,
        0x1000000008540001ULL,
        0xffffffff76c0001ULL,
        0xffffffff7680001ULL,
        0x1000000008e40001ULL,
        0xffffffff6fc0001ULL,
        0x1000000009240001ULL,
        0x1000000009300001ULL,
        0xffffffff6880001ULL,
        0x10000000098c0001ULL,
        0x1000000009c00001ULL,
        0xffffffff6340001ULL,
        0xffffffff5d40001ULL,
        0x100000000a8c0001ULL,
        0xffffffff54c0001ULL,
        0x100000000ab80001ULL,
        0x100000000aec0001ULL,
        0x100000000af40001ULL,
        0xffffffff4d40001ULL,
        0x100000000b300001ULL,
        0x100000000b400001ULL,
        0x100000000b7c0001ULL,
        0x100000000ba00001ULL,
        0xffffffff4380001ULL,
        0xffffffff3e80001ULL,
        0xffffffff37c0001ULL,
        0xffffffff36c0001ULL,
        0x100000000d400001ULL,
        0x100000000d4c0001ULL,
        0x100000000db00001ULL,
        0x100000000de00001ULL,
        0xffffffff2100001ULL,
        0x100000000df40001ULL,
        0x100000000e1c0001ULL,
        0xffffffff1d80001ULL,
        0xffffffff1cc0001ULL,
        0x100000000e3c0001ULL,
        0x100000000e480001ULL,
        0xffffffff1900001ULL,
        0xffffffff1740001ULL,
        0xffffffff15c0001ULL,
        0x100000000f140001ULL,
        0xffffffff0e80001ULL,
        0x100000000f200001ULL,
        0x100000000f3c0001ULL,
        0xfffffffeff80001ULL,
        0xfffffffeff40001ULL,
        0x1000000010100001ULL,
        0x1000000010a00001ULL,
        0xfffffffeefc0001ULL,
        0x1000000011340001ULL,
        0xfffffffee8c0001ULL
    };
    const int L = (2*LOGQ + 17 + 59) / 60;
	static CONV<L, LOGQ, N> *theConv = nullptr;
    static bool theConv_initialized = false;
    if(theConv_initialized == false) {
        theConv = new CONV<L, LOGQ, N>(Q_primes);
        theConv_initialized = true;
    }
    theConv->conv(A, B, C);
}