
#include "test.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  double Vr[N/2], Vi[N/2];
  double pt[N], pt_v_raw[N], pt_v[N], e[N], res[N];

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);
  double expected = N * square_sum_pt(pt) / 12.0;
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 100; ++r) {
  set_test_message(Vr, Vi);
  encode_raw(Vr, Vi, Delta, pt_v_raw);
  encode(Vr, Vi, Delta, pt_v);
  sub_pt(pt_v_raw, pt_v, e);
  conv(e, pt, res);

  double measured = square_sum_pt(res);
  std::cout << measured << std::endl;
  //std::cout << "Measured : " << measured << std::endl;
  //std::cout << "Expected : " << expected << std::endl;
  //std::cout << "Measured / Expected : " << (measured / expected) << std::endl;
  }
}
