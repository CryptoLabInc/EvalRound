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