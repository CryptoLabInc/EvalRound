
#include "define.h"
#include "message.h"
#include "plaintext.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  double Vr[N/2], Vi[N/2];
  int64_t pt[N];
  double e[N], res[N];


  set_test_message(zr, zi);
  set_test_message(Vr, Vi);
  encode(zr, zi, Delta, pt);
  encode_error(Vr, Vi, Delta, e);
  conv(e, pt, res);

  double measured = square_sum(res);
  double expected = N * square_sum(pt) / 12.0;
  std::cout << "Measured : " << measured << std::endl;
  std::cout << "Expected : " << expected << std::endl;
  std::cout << "Measured / Expected : " << (measured / expected) << std::endl; // currently 2
}
