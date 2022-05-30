#pragma once

#include <cstdint>


#ifdef _CONSOLE
extern "C"
{
    uint64_t mod_asm(uint64_t al, uint64_t ah, uint64_t q);
    uint64_t mul_lo_asm(uint64_t a, uint64_t b);
    uint64_t mul_hi_asm(uint64_t a, uint64_t b);
};

uint64_t mod(uint64_t al, uint64_t ah, uint64_t q) {
    return mod_asm(al, ah, q);
}

void mul_(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
    lo = mul_lo_asm(a, b);
    hi = mul_hi_asm(a, b);
}

#else

using u128 = unsigned __int128;

// compute (ah * beta + al) % q
uint64_t mod(uint64_t al, uint64_t ah, uint64_t q){
  u128 a = ah;
  u128 b = a << 64;
  u128 c = b + (u128) al;
  u128 d = c % (u128) q;
  return (uint64_t) d;
}

void mul_(uint64_t a, uint64_t b, uint64_t &lo, uint64_t &hi) {
  u128 a128 = a;
  u128 b128 = b;
  u128 c = a128 * b128;
  lo = c;
  hi = c >> 64;
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
  return add_mod(a, q - (b%q), q);
}

// compute (a * b) % q
uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t q) {
  uint64_t tempL, tempH;
  mul_(a, b, tempL, tempH);
  return mod(tempL, tempH, q);
}

// compute (a^n) % q
uint64_t pow_mod(uint64_t a, uint64_t n, uint64_t q) {
  if(n == 0)
    return 1;
  uint64_t res = pow_mod(a, n/2, q); // a^([n/2])
  res = mul_mod(res, res, q); // a^(2 * [n/2])
  if(n % 2== 1)
    res = mul_mod(a, res, q);
  return res;
}

// compute a^-1 (mod q) where q is prime
uint64_t inv_mod(uint64_t a, uint64_t q) {
  return pow_mod(a, q-2, q);
}