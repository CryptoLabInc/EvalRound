#include "experiment/rns_target.h"

int main()
{
  // ||Uz - \widetilde{Uz} || ~= sqrt(KN / 12) * (1 / Delta_A) * ||z|| * sqrt(D) * ||U_i||^{D-1}
  Message<LOGN> z, Uz, Uz_tilde, e;
  SparseDiagonal<N/2, 3> U0r[LOGN-1];
  SparseDiagonal<N/2, 3> U0i[LOGN-1];
  R_Q<LOGQ, N> pt, pt_Uz;

  const int D = LOGN-1; // number of matrices

  set_test_rounded_message(z, Delta);

  double matrix_norm_except_one = std::pow(sqrt(2), D-1);
  double expected = sqrt(3 * N / 12.0) * norm(z) / double (Delta) * sqrt(D) * matrix_norm_except_one;
  double bound = sqrt(D) * expected;
  std::cout << expected << std::endl;
  std::cout << bound << std::endl;
  std::cout << std::endl;
  
  // get Uz
  set_test_U0_matrix<LOGN>(U0r, U0i);
  matrix_vector_product<LOGN>(z, U0r[0], U0i[0], Uz);
  for(int i = 1; i < D; ++i) {
    Message<LOGN> temp(Uz);
    matrix_vector_product<LOGN>(temp, U0r[i], U0i[i], Uz);
  }

  // get Uz_tilde
  encode(z, Delta, pt);
  matrix_vector_product<LOGQ, LOGN>(pt, U0r[0], U0i[0], Delta, pt_Uz);
  for(int i = 1; i < D; ++i) {
    R_Q<LOGQ, N> temp(pt_Uz);
    matrix_vector_product<LOGQ, LOGN>(temp, U0r[i], U0i[i], Delta, pt_Uz);
  }
  decode_log<LOGQ, LOGN>(pt_Uz, LOGDELTA * (D + 1), Uz_tilde);

  // measure e
  sub(Uz, Uz_tilde, e);
  double measured = norm(e);
  std::cout << measured << std::endl;
}
