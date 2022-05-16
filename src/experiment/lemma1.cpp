
#include "experiment/simple.h"

#include <iostream>

int main()
{
  // ||e * pt|| ~= sqrt(N/12) * ||pt||
  Message <LOGN> z;
  SimplePlaintext<LOGN> pt;

  set_test_message(z);
  encode(z, Delta, pt);
  double expected = N /12.0 * square_sum(pt);
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 100; ++r) {
    Message <LOGN> V;
    SimplePlaintext<LOGN> pt_v_raw, pt_v, e, res;

    set_random_message(V);
    encode_raw(V, Delta, pt_v_raw);
    encode(V, Delta, pt_v);
    sub(pt_v_raw, pt_v, e);
    conv(e, pt, res);

    double measured = square_sum(res);
    std::cout << measured << std::endl;
  }
}
