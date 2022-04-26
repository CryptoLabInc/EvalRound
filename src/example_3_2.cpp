#include "common.h"
#include "HEAAN/HEAAN.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  R_Q<LOGQ, N> pt;
  R_Q<LOGQ_UP, N> pt_up;
  int s[N];
  R_Q_square<LOGQ, N> ct;
  R_Q_square<LOGQ_UP, N> ct_up;

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);
  decode(pt, Delta, zr, zi);
  HEAAN<LOGQ, N>::keygen(h, s);
  HEAAN<LOGQ, N>::enc(pt, s, ct);
  mod_raise<LOGQ, LOGQ_UP, N> (pt, pt_up);
  HEAAN<LOGQ_UP, N>::dec(ct_up, s, pt_up);

  // now what? help!

  // check distribution of pt+qI
  // gen normal distribution
  // check similarity
}
