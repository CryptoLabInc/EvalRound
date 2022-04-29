
#include "test.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  double Ar[K][N/2], Ai[K][N/2];
  double Azr[N/2], Azi[N/2];

  double pt[N];
  double pt_Az_raw[N], pt_Az[N];
  double e_Az[N];

  // get pt_Az_raw
  set_test_rounded_message(zr, zi);
  set_test_matrix(Ar, Ai);
  matrix_vector_product(zr, zi, Ar, Ai, Azr, Azi);
  encode(Azr, Azi, DeltaSq, pt_Az_raw);

  // get pt_Az
  encode(zr, zi, Delta, pt);
  matrix_vector_product(pt, Ar, Ai, pt_Az);

  // get e_Az
  sub_pt(pt_Az_raw, pt_Az, e_Az);

  // sanity check
  double er[N/2], ei[N/2];
  decode(e_Az, DeltaSq, er, ei);
  print_pt("er (for sanity check)", er);

  double measured = square_sum_pt(e_Az);
  double expected = K * N * square_sum_pt(pt) / 12.0;
  double bound = K * expected;
  std::cout << "Measured : " << measured << std::endl;
  std::cout << "Expected : " << expected << std::endl;
  std::cout << "Bound : " << bound << std::endl;
  std::cout << "Measured / Expected : " << (measured / expected) << std::endl;
}
