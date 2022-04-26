
#include "HEAAN/HEAAN.h"
#include "common.h"

#include <iostream>

int main()
{
  double zr[N/2], zi[N/2];
  SparseDiagonal<N/2, K> Ar, Ai;
  double Azr[N/2], Azi[N/2];

  R_Q<LOGQ, N> pt;
  R_Q<LOGQ, N> pt_v[K];

  set_test_message(zr, zi);
  set_test_matrix(Ar, Ai);
  R_Q<LOGQ, N> pt_Az;
  matrix_vector_product(zr, zi, Ar, Ai, Azr, Azi);

  encode(zr, zi, Delta, pt);
  for(int k = 0; k < K; ++k)
    encode(Ar.vec[k], Ai.vec[k], Delta, pt_v[k]);
  encode(Azr, Azi, DeltaSq, pt_Az);

  R_Q<LOGQ, N> pt_conv;
  conv<LOGQ, N>(pt_v[0], pt, pt_conv);

  std::cout << pt_Az.square_sum() << std::endl;
  pt_Az.print_signed();
  std::cout << pt_conv.square_sum() << std::endl;
  pt_conv.print_signed();
 
  // pt_A * pt_z
/*  R_Q<LOGQ, N> pt_Az_tilde;
  double z_rot_r[N/2];
  double z_rot_i[N/2];
  R_Q<LOGQ, N> pt_rot;
  R_Q<LOGQ, N> pt_conv;
  pt_Az_tilde.setzero();
  for(int k = 0; k < K; ++k) {
    rotate_message(zr, z_rot_r, k);
    rotate_message(zi, z_rot_i, k);
    encode(z_rot_r, z_rot_i, Delta, pt_rot);
    conv<LOGQ, N>(pt_v[k], pt_rot, pt_conv);
    pt_Az_tilde += pt_conv;
  }

  R_Q<LOGQ, N> e_Az = pt_Az;
  e_Az -= pt_Az_tilde;
  uint64_t e_Az_sq_sum = e_Az.square_sum(); 
  uint64_t pt_sq_sum = pt.square_sum();

  std::cout << pt_Az.square_sum() << std::endl;
  std::cout << pt_Az_tilde.square_sum() << std::endl;
  std::cout << e_Az_sq_sum << " " << pt_sq_sum << std::endl;
*/
  // measure ptxt side rounding error
  // measure pt norm
}
