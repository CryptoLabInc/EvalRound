
#include "define.h"
#include "message.h"
#include "plaintext.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  double Ar[K][N/2], Ai[K][N/2];
  double Azr[N/2], Azi[N/2];

  int64_t pt[N];
  int64_t pt_v[K][N];
  int64_t pt_Az[N];
  int64_t pt_Az_tilde[N];
  int64_t e_Az[N];

  set_test_message(zr, zi);
  set_test_matrix(Ar, Ai);
  matrix_vector_product(zr, zi, Ar, Ai, Azr, Azi);
  encode(Azr, Azi, DeltaSq, pt_Az);

  encode(zr, zi, Delta, pt);
  matrix_vector_product(pt, Ar, Ai, pt_Az_tilde);

  double Azr_tilde[N/2], Azi_tilde[N/2];
  decode(pt_Az_tilde, DeltaSq, Azr_tilde, Azi_tilde);

  sub(pt_Az, pt_Az_tilde, e_Az);
  double er[N/2], ei[N/2];
  decode(e_Az, DeltaSq, er, ei);
  print("er (for sanity check)", er);

  double measured = square_sum(e_Az);
  double expected = K * N * square_sum(pt) / 12.0;
  double bound = K * expected;
  std::cout << "Measured : " << measured << std::endl;
  std::cout << "Expected : " << expected << std::endl;
  std::cout << "Bound : " << bound << std::endl;
  std::cout << "Measured / Expected : " << (measured / expected) << std::endl; // currently 2
}
