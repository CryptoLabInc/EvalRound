#include "test.h"

#include <cmath>

void example1() {
  double zr[N/2], zi[N/2];
  double Ar[K][N/2], Ai[K][N/2];
  double Br[K][N/2], Bi[K][N/2];
  double Azr[N/2], Azi[N/2];
  double BAzr[N/2], BAzi[N/2];
  double BAzr_tilde[N/2], BAzi_tilde[N/2];
  double er[N/2], ei[N/2];

  double pt[N];
  double pt_Az[N];
  double pt_BAz[N];

  // get Az
  set_test_rounded_message(zr, zi);
  set_test_matrix(Ar, Ai);
  set_test_matrix(Br, Bi);
  matrix_vector_product(zr, zi, Ar, Ai, Azr, Azi);
  matrix_vector_product(Azr, Azi, Br, Bi, BAzr, BAzi);

  // get Az_tilde
  encode(zr, zi, Delta, pt);
  matrix_vector_product(pt, Ar, Ai, pt_Az);
  matrix_vector_product(pt_Az, Br, Bi, pt_BAz);
  decode(pt_BAz, DeltaTr, BAzr_tilde, BAzi_tilde);

  // get e
  sub(BAzr, BAzr_tilde, er);
  sub(BAzi, BAzi_tilde, ei);

  // sanity check
  print("er (for sanity check)", er);
}

int main()
{
  example1();
}
