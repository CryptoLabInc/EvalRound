
#include "define.h"
#include "message.h"
#include "plaintext.h"

#include <iostream>
#include <cmath>

int main()
{
  double zr[N/2], zi[N/2];
  int64_t pt[N];

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);

  double measured = norm(pt) / norm(zr, zi);
  double expected = (double) Delta / sqrt(N/2);

  print("zr", zr);
  print("pt", pt);

  std::cout << "Measured : " << measured << std::endl;
  std::cout << "Expected : " << expected << std::endl;



  // measure decoded diff from 3.3
  // compare
}
