#pragma once

#include "crt.h"
#include "ntt.h"
#include "mod.h"

#include <memory>

template<int L, int N>
struct CONV{
    CRT<L> crt;
    std::unique_ptr<NTT<N>> ntt[L];

    CONV(const uint64_t q[L]);

    template<int LOGQ>
    void conv( const R_Q<LOGQ,N>& A,
	       const R_Q<LOGQ,N>& B, R_Q<LOGQ,N>& C );

    private:
        template<int LOGQ>
        void icrt_polynomial(const R_Q<LOGQ, N>& A, uint64_t A_rns[L][N]);

        template<int LOGQ>
        void crt_polynomial(const uint64_t A_rns[L][N], R_Q<LOGQ, N> &A);
};

template<int L, int N>
CONV<L, N>::CONV(const uint64_t q[L]) : crt{q} {
    for(int i = 0; i < L; ++i) {
        ntt[i] = std::make_unique<NTT<N>>(q[i]);
    }
}

template<int L, int N>
template<int LOGQ>
void CONV<L, N>::conv( const R_Q<LOGQ,N>& A,
	       const R_Q<LOGQ,N>& B, R_Q<LOGQ,N>& C ){
    uint64_t A_rns[L][N], B_rns[L][N], C_rns[L][N];

    icrt_polynomial(A, A_rns);

    for(int i = 0; i < L; ++i) {
        ntt[i]->ntt(A_rns[i]);
        ntt[i]->ntt(B_rns[i]);
        for(int j = 0; j < N; ++j) {
            C_rns[i][j] = mul_mod(A_rns[i][j], B_rns[i][j], ntt[i]->q);
        }
        ntt[i]->intt(C_rns[i]);
    }

    crt_polynomial(C_rns, C);
}

template<int L, int N>
template<int LOGQ>
void CONV<L, N>::icrt_polynomial(const R_Q<LOGQ, N>& A, uint64_t A_rns[L][N]) {
    for(int i = 0; i < N; ++i) {
        Z_Q<LOGQ> Ai(A[i]);
        uint64_t Ai_arr[L], Ai_rns[L];

        for(int j = 0; j < L; ++j){
            Ai_arr[j] = j < ((63 + LOGQ) / 64) ? Ai[j] : 0;
        }
        crt.icrt(Ai_arr, Ai_rns);

        for(int j = 0; j < L; ++j){
            A_rns[j][i] = Ai_rns[j];
        }
    }
}

template<int L, int N>
template<int LOGQ>
void CONV<L, N>::crt_polynomial(const uint64_t A_rns[L][N], R_Q<LOGQ, N> &A) {
    for(int i = 0; i < N; ++i) {
        uint64_t Ai_rns[L], Ai[L];
        for(int j = 0; j < L; ++j) {
            Ai_rns[j] = A_rns[j][i];
        }
        crt.crt(Ai_rns, Ai);
        A[i].setzero();
        for(int j = 0; j < (63 + LOGQ) / 64; ++j) {
            A[i][j] = Ai[j];
        }
    }
}