#pragma once

#include "HEAAN/Z_Q.h"

#include <cstdint>


#ifdef _CONSOLE
#else

// compute (ah * beta + al) % q
uint64_t mod(uint64_t al, uint64_t ah, uint64_t q){
  using u128 = unsigned __int128;
  u128 a = ah;
  u128 b = a << 64;
  u128 c = b + (u128) al;
  u128 d = c % (u128) q;
  return (uint64_t) d;
}
#endif

// compute (a + b) % q
uint64_t add_mod(uint64_t a, uint64_t b, uint64_t q) {
  uint64_t tempL, tempH;
  tempL = a + b;
  tempH = tempL < a;
  return mod(tempL, tempH, q);
}

// compute (a - b) % q
uint64_t sub_mod(uint64_t a, uint64_t b, uint64_t q) {
  return add_mod(a, q - b, q);
}

// compute (a * b) % q
uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t q) {
  uint64_t tempL, tempH;
  mul(a, b, tempL, tempH);
  return mod(tempL, tempH, q);
}

// compute (a^n) % q
uint64_t power_mod(uint64_t a, uint64_t n, uint64_t q) {
  if(n == 0)
    return 1;
  uint64_t res = power_mod(a, n/2, q); // a^([n/2])
  res = mul_mod(res, res, q); // a^(2 * [n/2])
  if(n % 2== 1)
    res = mul_mod(a, res, q);
  return res;
}

// compute a^-1 (mod q) where q is prime
uint64_t inv_mod(uint64_t a, uint64_t q) {
  return power_mod(a, q-1, q);
}