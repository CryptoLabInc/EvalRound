
#include "experiment/rns_target.h"

#include <iostream>

int main()
{
  // ||Az - \widetilde{Az} || ~= sqrt(DiagN / 12) * (1 / Delta_A) * ||z||
  Message<LOGN> z;
  R_Q<LOGQ, N> pt;

  const int Diag = 16; // number of diagonals on matrix

  set_test_rounded_message(z, Delta);
  encode(z, Delta, pt);

  double expected = sqrt(Diag * N / 12.0) * norm(z) / (double) Delta;
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 100; ++r) {
    Message<LOGN> A[Diag];
    Message<LOGN> Az, Az_tilde, e;
    R_Q<LOGQ, N> pt_Az;

    // get Az
    set_random_matrix<LOGN, Diag>(A);
    matrix_vector_product<LOGN, Diag>(z, A, Az);

    // get Az_tilde
    matrix_vector_product<LOGQ, LOGN, Diag>(pt, A, Delta, pt_Az);
    decode_log(pt_Az, 2*LOGDELTA, Az_tilde);

    // measure e
    sub(Az, Az_tilde, e);
    double measured = norm(e);
    std::cout << measured << std::endl;
  }
}
