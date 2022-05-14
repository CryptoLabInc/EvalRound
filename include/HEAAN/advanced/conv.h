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
    // 20 61bit primes, mod(q[i], 2^17) = 1
    static constexpr uint64_t Q_primes[] = {
	    0x1fffffffffe00001ULL,
	    0x1fffffffffc80001ULL,
	    0x2000000000460001ULL,
	    0x1fffffffffb40001ULL,
	    0x2000000000500001ULL,
        0x2000000000620001ULL,
        0x20000000007c0001ULL,
        0x2000000000a00001ULL,
        0x2000000000b00001ULL,
        0x1fffffffff500001ULL,
        0x1fffffffff420001ULL,
        0x1fffffffff380001ULL,
        0x2000000000da0001ULL,
        0x2000000000ec0001ULL,
        0x2000000000f80001ULL,
        0x1fffffffff000001ULL,
        0x1ffffffffef00001ULL,
        0x2000000001120001ULL,
        0x1ffffffffee80001ULL,
        0x2000000001480001ULL,
        0x1ffffffffeb40001ULL,
        0x1ffffffffe780001ULL,
        0x2000000001960001ULL,
        0x1ffffffffe600001ULL,
        0x1ffffffffe4c0001ULL,
        0x2000000001be0001ULL,
        0x2000000001dc0001ULL,
        0x2000000001e40001ULL,
        0x2000000002060001ULL,
        0x1ffffffffdf40001ULL,
        0x1ffffffffdce0001ULL,
        0x20000000023e0001ULL,
        0x1ffffffffdb20001ULL,
        0x1ffffffffdac0001ULL,
        0x1ffffffffda40001ULL,
        0x2000000002660001ULL,
        0x2000000002680001ULL,
        0x20000000026c0001ULL,
        0x2000000002740001ULL,
        0x2000000002800001ULL,
        0x1ffffffffd7a0001ULL,
        0x2000000002c20001ULL,
        0x2000000002ec0001ULL,
        0x2000000003020001ULL,
        0x20000000033e0001ULL,
        0x2000000003500001ULL,
        0x20000000036a0001ULL,
        0x2000000003940001ULL,
        0x1ffffffffc680001ULL,
        0x2000000003a60001ULL,
        0x1ffffffffc000001ULL,
        0x20000000041e0001ULL,
        0x2000000004240001ULL,
        0x2000000004280001ULL,
        0x1ffffffffb880001ULL,
        0x1ffffffffb7c0001ULL,
        0x20000000048e0001ULL,
        0x2000000004ac0001ULL,
        0x2000000004b40001ULL,
        0x2000000004ba0001ULL,
        0x1ffffffffb300001ULL,
        0x1ffffffffb1e0001ULL,
        0x1ffffffffb1c0001ULL,
        0x1ffffffffb0a0001ULL,
        0x2000000005080001ULL,
        0x1ffffffffaf20001ULL,
        0x1ffffffffadc0001ULL
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