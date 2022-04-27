#include "define.h"
#include "message.h"
#include "plaintext.h"

#include <iostream>
#include <cmath>

void rotate_test(){
  double zr[N/2], zi[N/2];
  double zr_rot1[N/2], zi_rot1[N/2];
  double zr_rot2[N/2], zi_rot2[N/2];
  set_test_message(zr, zi);
  print("zr", zr);

  int64_t pt[N], pt_rot[N];
  encode(zr, zi, Delta, pt);
  for(int i = 0; i < 10; ++i) {
    std::cout << i << std::endl;
    rotate(zr, zr_rot1, i);
    print("zr_rot1", zr_rot1);

    rotate(pt, pt_rot, i);
    decode(pt_rot, Delta, zr_rot2, zi_rot2);
    print("zr_rot2", zr_rot2);
  }
}

int main()
{
  double zr[N/2], zi[N/2];
  double Ar[K][N/2], Ai[K][N/2];
  double Br[K][N/2], Bi[K][N/2];
  double Azr[N/2], Azi[N/2];
  double BAzr[N/2], BAzi[N/2];

  set_test_message(zr, zi);
  set_test_matrix(Ar, Ai);
  set_test_matrix(Br, Bi);
  matrix_vector_product(zr, zi, Ar, Ai, Azr, Azi);
  matrix_vector_product(Azr, Azi, Br, Bi, BAzr, BAzi);

  int64_t pt[N];
  int64_t pt_Az[N];
  int64_t pt_BAz[N];
  double BAzr_tilde[N/2], BAzi_tilde[N/2];
  double e_BAzr[N], e_BAzi[N];

  encode(zr, zi, Delta, pt);
  matrix_vector_product(pt, Ar, Ai, pt_Az);
  matrix_vector_product(pt_Az, Br, Bi, pt_BAz);
  decode(pt_BAz, DeltaTr, BAzr_tilde, BAzi_tilde);

  sub(BAzr, BAzr_tilde, e_BAzr);
  sub(BAzi, BAzi_tilde, e_BAzi);
  print("BAzr", BAzr);
  print("BAzr_tilde", BAzr_tilde);
  print("e_BAzr (sanity check)", e_BAzr);
}
