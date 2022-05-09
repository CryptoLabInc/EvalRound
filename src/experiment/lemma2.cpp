
#include "experiment/simple.h"

#include <iostream>

int main()
{
  // || \sigma_k e * pt || ~= sqrt( K N /  12) * || pt||
  Message <LOGN> z;
  SimplePlaintext<LOGN> pt;

  set_test_message(z);
  encode(z, Delta, pt);
  double expected = sqrt(K * N / 12.0) * norm(pt);
  std::cout << expected << std::endl;
  std::cout << std::endl;
  
  for(int r = 0; r < 100; ++r) {
    Message <LOGN> V;
    SimplePlaintext<LOGN> pt_v_raw, pt_v, pt_rot, pt_conv, e, res;
    for(int k = 0; k < K; ++k) {
      set_random_message(V);
      encode_raw(V, Delta, pt_v_raw);
      encode(V, Delta, pt_v);
      sub(pt_v_raw, pt_v, e);
      rotate(pt, pt_rot, k);
      conv(e, pt, pt_conv);
      add(pt_conv, res, res);
    }

    double measured = norm(res);
    std::cout << measured << std::endl;
  }
}
