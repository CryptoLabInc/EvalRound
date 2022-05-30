#include "experiment/section3.h"

#include <iostream>

template<int LOGN>
void measure_conv_error()
{
  std::cout << "Measuring error on convolution" << std::endl;
  std::cout << "LOGN : " << LOGN << std::endl;
  const int N = 1 << LOGN;

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

int main()
{
  measure_conv_error<15>();
  measure_conv_error<16>();
}