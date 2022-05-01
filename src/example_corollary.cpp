#include "SIMPLE/test.h"

#include <cmath>

void example() {
  // multi matrix version of theorem1
  double zr[N/2], zi[N/2];
  SparseDiagonal<N/2, 3> U0r[LOGN-1];
  SparseDiagonal<N/2, 3> U0i[LOGN-1];
  double Uzr[D][N/2], Uzi[D][N/2];

  // get Uz
  set_test_rounded_message(zr, zi);
  set_test_U0_matrix(U0r, U0i);
  matrix_vector_product(zr, zi, U0r[0], U0i[0], Uzr[0], Uzi[0]);
  for(int i = 1; i < D; ++i) {
    matrix_vector_product(Uzr[i-1], Uzi[i-1], U0r[i], U0i[i], Uzr[i], Uzi[i]);
  }

  double matrix_norm_except_one = std::pow(sqrt(2), D-1);
  double expected = sqrt(3 * N / 12.0) * norm(zr, zi) / double (Delta) * D * matrix_norm_except_one;
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 1; ++r) {
    double Uzr_tilde[N/2], Uzi_tilde[N/2];
    double er[N/2], ei[N/2];
    double pt[N], pt_Uz[D][N];
    encode(zr, zi, Delta, pt);

    // get Uz_tilde
    matrix_vector_product(pt, U0r[0], U0i[0], pt_Uz[0]);
    for(int i = 1; i < D; ++i) {
      matrix_vector_product(pt_Uz[i-1], U0r[i], U0i[i], pt_Uz[i]);
    }
    decode(pt_Uz[D-1], DeltaTotal, Uzr_tilde, Uzi_tilde);

    // get e
    sub(Uzr[D-1], Uzr_tilde, er);
    sub(Uzi[D-1], Uzi_tilde, ei);

    // sanity check
    //print("er (for sanity check)", er);

    double measured = norm(er, ei);
    std::cout << measured << std::endl;
    //double bound = sqrt(3) * expected;
    //std::cout << "Measured : " << measured << std::endl;
    //std::cout << "Expected : " << expected << std::endl;
    //std::cout << "Bound : " << bound << std::endl;
    //std::cout << "Measured / Expected : " << (measured / expected) << std::endl;
  }
}

void matrix_norm_test() {
  double zr[N/2], zi[N/2];
  double Azr[N/2], Azi[N/2];
  SparseDiagonal<N/2, 3> U0r[LOGN-1];
  SparseDiagonal<N/2, 3> U0i[LOGN-1];

  set_test_message(zr, zi);
  set_test_U0_matrix(U0r, U0i);
  for(int i = 0; i < LOGN-1; ++i) {
    matrix_vector_product(zr, zi, U0r[i], U0i[i], Azr, Azi);
    double Az_sq = square_sum(Azr, Azi);
    double z_sq = square_sum(zr, zi);
    std::cout << (Az_sq / z_sq) << std::endl;
  }
}

int main()
{
  example();
}
