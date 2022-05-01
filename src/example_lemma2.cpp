
#include "test.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  double Vr[N/2], Vi[N/2];
  double pt[N], pt_v_raw[N], pt_v[N], pt_rot[N], pt_conv[N], e[N], res[N];

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);
  double expected = sqrt(K * N / 12.0) * norm_pt(pt);
  std::cout << expected << std::endl;
  std::cout << std::endl;
  
  for(int r = 0; r < 100; ++r) {
    set_zero(res);
    for(int k = 0; k < K; ++k) {
      set_random_message(Vr, Vi);
      encode_raw(Vr, Vi, Delta, pt_v_raw);
      encode(Vr, Vi, Delta, pt_v);
      sub_pt(pt_v_raw, pt_v, e);
      rotate(pt, pt_rot, k);
      conv(e, pt, pt_conv);
      add(pt_conv, res, res);
    }

    double measured = norm_pt(res);
    std::cout << measured << std::endl;
    //std::cout << "Measured : " << measured << std::endl;
    //std::cout << "Expected : " << expected << std::endl;
    //std::cout << "Measured / Expected : " << (measured / expected) << std::endl;
  }
}
