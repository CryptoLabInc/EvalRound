
#include "SIMPLE/test.h"

#include <iostream>

int main()
{
  // ||e * pt|| ~= sqrt(N/12) * ||pt||
  double zr[N/2], zi[N/2];
  double Vr[N/2], Vi[N/2];
  double pt[N], pt_v_raw[N], pt_v[N], e[N], res[N];

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);
  double expected = sqrt(N /12.0) * norm_pt(pt);
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 100; ++r) {
    set_random_message(Vr, Vi);
    encode_raw(Vr, Vi, Delta, pt_v_raw);
    encode(Vr, Vi, Delta, pt_v);
    sub_pt(pt_v_raw, pt_v, e);
    conv(e, pt, res);

    double measured = norm_pt(res);
    std::cout << measured << std::endl;
    //std::cout << "Measured : " << measured << std::endl;
    //std::cout << "Expected : " << expected << std::endl;
    //std::cout << "Measured / Expected : " << (measured / expected) << std::endl;
  }
}
