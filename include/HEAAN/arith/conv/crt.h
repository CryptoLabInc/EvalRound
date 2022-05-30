#pragma once

#include "HEAAN/arith/Z_Q.h"
#include "mod.h"

namespace {
  template<int LOGQ>
  bool cmpge(const Z_Q <LOGQ> &A, const Z_Q <LOGQ> &B){
    for (int i = (63 + LOGQ) / 64 - 1; i >= 0; --i) {
      if(A[i] > B[i]) return true;
      else if(A[i] < B[i]) return false;
    }
    return true;
  }
}

template<int L, int LOGQ>
struct CRT {
  Z_Q<64*(L+1)> Q; // Original modulus
  uint64_t q[L]; // Decomposed modulus s.t. Q = \Pi q[i]

  // if a = a_rns[i] (mod q[i]),
  // a = crt(a_rns) = \Sigma a_rns[i] (Q/q[i])^-1 * (Q/q[i]) (mod Q)
  Z_Q<64*(L+1)> Q_over_q[L]; // Q/q[i]
  uint64_t inv_mod_Q_over_q[L]; // (Q/q[i])^-1 (mod q[i])

  CRT(const uint64_t q[L]);
  void crt(const uint64_t a_rns[L], Z_Q<LOGQ> &A) const;
  void icrt(const Z_Q<LOGQ> &A, uint64_t a_rns[L]) const;
};

template<int L, int LOGQ>
CRT<L, LOGQ>::CRT(const uint64_t q[L]) {
  for(int i = 0; i < L; ++i)
    this->q[i] = q[i];
  
  // compute Q
  Q.setzero();
  Q[0] = 1;
  for(int i = 0; i < L; ++i)
    Q *= q[i];

  // compute Q_over_q
  for(int i = 0; i < L; ++i){
    Q_over_q[i].setzero();
    Q_over_q[i][0] = 1;
    for(int j = 0; j < L; ++j){
      if(i != j)
        Q_over_q[i] *= q[j];
    }
  }
  
  // compute inv_mod_Q
  for(int i = 0; i < L; ++i) {
    uint64_t temp = 1;
    for(int j = 0; j < L; ++j) {
      if(i != j)
        temp = mul_mod(temp, q[j], q[i]);
    }
    inv_mod_Q_over_q[i] = inv_mod(temp, q[i]);
  }
}


template<int L, int LOGQ>
void CRT<L, LOGQ>::crt(const uint64_t a_rns[L], Z_Q<LOGQ> &A) const {
    Z_Q<64*(L+1)> sum, neg_sum;
    sum.setzero();
    for(int i =0; i < L; ++i) {
      // coeff = (a_rns[i] * (Q/q[i])^-1) % q[i]
      uint64_t coeff = mul_mod(a_rns[i], inv_mod_Q_over_q[i], q[i]);
      // sum = \Sigma coeff * Q_over_q[i];
      Z_Q<64*(L+1)> temp(Q_over_q[i]);
      temp *= coeff;
      sum += temp;
    }
    // sum = sum % Q
    while(cmpge(sum, Q)) {
      sum -= Q;
    }

    // handling sign
    neg_sum = Q;
    neg_sum -= sum;
    bool is_negative = cmpge(sum, neg_sum);
    if(is_negative) {
      resize(neg_sum, A);
      A.negate();
    } else {
      resize(sum, A);
    }
  }

template<int L, int LOGQ>
void CRT<L, LOGQ>::icrt(const Z_Q<LOGQ> &A, uint64_t a_rns[L]) const {
  bool is_negative = A.is_bigger_than_halfQ();
  if(is_negative) {
    Z_Q<LOGQ> A_abs(A);
    A_abs.negate();

    int len = A_abs.get_length();
    for(int i = 0; i< L; ++i){
      // compute a % q[i] where a = \Sigma a[i]*beta^i
      uint64_t ah = A_abs[len-1];
      for(int j = len-2; j>=0; --j) {
        ah  = mod(A_abs[j], ah, q[i]);
      }
      a_rns[i] = q[i] - ah;
    }
  } else {
    int len = A.get_length();
    for(int i = 0; i< L; ++i){
      // compute a % q[i] where a = \Sigma a[i]*beta^i
      uint64_t ah = A[len-1];
      for(int j = len-2; j>=0; --j) {
        ah  = mod(A[j], ah, q[i]);
      }
      a_rns[i] = ah;
    }
  }
}