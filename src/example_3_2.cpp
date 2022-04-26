#include "common.h"
#include "HEAAN/HEAAN.h"

#include <vector>

int main()
{
  double *zr = new double(N/2), *zi = new double(N/2);
  R_Q<LOGQ, 1> pt;

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);

  // set input message
  // encode it into pt
  // randomly gen sk
  // encrypt
  // mod raise
  // check distribution of pt+qI
  // gen normal distribution
  // check similarity

  delete zr, zi;
}
