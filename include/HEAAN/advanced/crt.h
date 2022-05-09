#pragma once

#include "HEAAN/Z_Q.h"
#include "mod.h"

#include <iostream>

template<int L>
struct CRT {
  Z_Q<64*(L+1)> Q; // Original modulus
  uint64_t q[L]; // Decomposed modulus s.t. Q = \Pi q[i]

  // if a = a_rns[i] (mod q[i]),
  // a = crt(a_rns) = \Sigma a_rns[i] (Q/q[i])^-1 * (Q/q[i]) (mod Q)
  Z_Q<64*(L+1)> Q_over_q[L]; // Q/q[i]
  uint64_t inv_mod_Q_over_q[L]; // (Q/q[i])^-1 (mod q[i])

  CRT(const uint64_t q[L],
  const uint64_t Q[L],
  const uint64_t inv_mod_Q_over_q[L],
  const uint64_t Q_over_q[L][L]) {
    for(int i = 0; i < L; ++i)
      this->q[i] = q[i];
    for(int i = 0; i < L; ++i)
      this->Q[i] = Q[i];
    this->Q[L] = 0;
    for(int i = 0; i < L; ++i)
      this->inv_mod_Q_over_q[i] = inv_mod_Q_over_q[i];
    for(int i = 0; i < L; ++i){
      for(int j = 0; j < L; ++j)
        this->Q_over_q[i][j] = Q_over_q[i][j];
      this->Q_over_q[i][L] = 0;
    }
  }

  void crt(const uint64_t a_rns[L], uint64_t a[L]);
  void icrt(const uint64_t a[L], uint64_t a_rns[L]);

  private:
    bool cmpge(const Z_Q <64*(L+1)> &A, const Z_Q <64*(L+1)> &B);
};


template<int L>
void CRT<L>::crt(const uint64_t a_rns[L], uint64_t a[L]) {
    Z_Q<64*(L+1)> sum;
    sum.setzero();
    Z_Q<64*(L+1)> temp;
    for(int i =0; i < L; ++i) {
      // coeff = (a_rns[i] * (Q/q[i])^-1) % q[i]
      uint64_t coeff = mul_mod(a_rns[i], inv_mod_Q_over_q[i], q[i]);

      // sum = \Sigma coeff * Q_over_q[i];
      temp = Q_over_q[i];
      temp *= coeff;
      sum += temp;
    }

    // sum = sum % Q
    while(cmpge(sum, Q))
      sum -= Q;

    // return
    for(int i = 0; i < L; ++i)
      a[i] = sum[i];
  }

template<int L>
void CRT<L>::icrt(const uint64_t a[L], uint64_t a_rns[L]){
  for(int i = 0; i<L; ++i){
    // compute a % q[i] where a = \Sigma a[i]*beta^i
    uint64_t ah = a[L-1];
    for(int j = L-2; j>=0; --j) {
      ah  = mod(a[j], ah, q[i]);
    }
    a_rns[i] = ah;
  }
}

template<int L>
bool CRT<L>::cmpge(const Z_Q <64*(L+1)> &A, const Z_Q <64*(L+1)> &B){
  for (int i =L; i >= 0; --i) {
    if(A[i] > B[i]) return true;
    else if(A[i] < B[i]) return false;
  }
  return true;
}
