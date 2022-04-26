#include "common.h"
#include "HEAAN/HEAAN.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  R_Q<LOGQ, N> pt;

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);
  decode(pt, Delta, zr, zi);

  // set input message
  // encode it into pt
  // randomly gen sk
  // encrypt
  // mod raise
  // check distribution of pt+qI
  // gen normal distribution
  // check similarity
}
