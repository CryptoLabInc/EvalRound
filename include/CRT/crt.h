#include <cstdint>
#include <iostream>

#include "HEAAN/Z_Q.h"

using u128 = unsigned __int128;

uint64_t upper(u128 a) {
  return (uint64_t) (a >> 64);
}

uint64_t lower(u128 a) {
  return (uint64_t) a;
}

void print(u128 a) {
  std::cout << upper(a) << " " << lower(a) << std::endl;
}

uint64_t mod(uint64_t ah, uint64_t al, uint64_t q){
  u128 a = ah;
  u128 b = a << 64;
  u128 c = b + (u128) al;
  u128 d = c % (u128) q;
  return (uint64_t) d;
}

template<int L>
struct CRT {
  uint64_t q[L];
  Z_Q<64*(L+1)> Q;
  uint64_t inv_mod_Q_over_q[L];
  Z_Q<64*(L+1)> Q_over_q[L];

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

  void crt(const uint64_t a_rns[L], uint64_t a[L]) {
    Z_Q<64*(L+1)> sum;
    sum.setzero();
    Z_Q<64*(L+1)> temp;
    for(int i =0; i < L; ++i) {
      uint64_t tempL, tempH;
      mul(a_rns[i], inv_mod_Q_over_q[i], tempL, tempH);
      temp = Q_over_q[i];
      temp *= mod(tempH, tempL, q[i]);
      sum += temp;
    }
    while(cmpge(sum, Q))
      sum -= Q;
    for(int i = 0; i < L; ++i)
      a[i] = sum[i];
  }

  void icrt(const uint64_t a[L], uint64_t a_rns[L]);

  private:
    bool cmpge(const Z_Q <64*(L+1)> &A, const Z_Q <64*(L+1)> &B);
};

template<int L>
void CRT<L>::icrt(const uint64_t a[L], uint64_t a_rns[L]){
  for(int i = 0; i<L; ++i){
    uint64_t ah = a[L-1];
    for(int j = L-2; j>=0; --j) {
      ah  = mod(ah, a[j], q[i]);
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
