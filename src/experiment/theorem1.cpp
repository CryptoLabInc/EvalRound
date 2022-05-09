
#include "experiment/big.h"

#include <iostream>

int main()
{
  // ||Az - \widetilde{Az} || ~= sqrt(KN / 12) * (1 / Delta_A) * ||z||
  Message<LOGN> z;
  R_Q<LOGQ, N> pt;

  set_test_rounded_message(z, Delta);
  encode(z, Delta, pt);

  double expected = sqrt(K * N / 12.0) * norm(z) / (double) Delta;
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 100; ++r) {
    Message<LOGN> A[K];
    Message<LOGN> Az, Az_tilde, e;
    R_Q<LOGQ, N> pt_Az;

    // get Az
    set_random_matrix<LOGN, K>(A);
    matrix_vector_product<LOGN, K>(z, A, Az);

    // get Az_tilde
    matrix_vector_product<LOGQ, LOGN, K>(pt, A, Delta, pt_Az);
    decode_log(pt_Az, 2*LOGDELTA, Az_tilde);

    // measure e
    sub(Az, Az_tilde, e);
    double measured = norm(e);
    std::cout << measured << std::endl;
  }
}
