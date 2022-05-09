#include <cstdint>


#ifdef _CONSOLE
#else
uint64_t mod(uint64_t ah, uint64_t al, uint64_t q){
  using u128 = unsigned __int128;
  u128 a = ah;
  u128 b = a << 64;
  u128 c = b + (u128) al;
  u128 d = c % (u128) q;
  return (uint64_t) d;
}
#endif