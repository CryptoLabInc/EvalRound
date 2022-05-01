
#include "SIMPLE/test.h"

#include <iostream>

void example1() {
  double zr[N/2], zi[N/2];
  double pt[N];

  set_test_message(zr, zi);
  encode(zr, zi, Delta, pt);

  double measured = norm_pt(pt) / norm(zr, zi);
  double expected = (double) Delta / sqrt(N/2);

  std::cout << "Measured : " << measured << std::endl;
  std::cout << "Expected : " << expected << std::endl;
}

int main()
{
  // ||Az - \widetilde{Az} || ~= sqrt(KN / 12) * (1 / Delta_A) * ||pt||
  double zr[N/2], zi[N/2];
  double pt[N];

  // get pt
  set_test_rounded_message(zr, zi);
  encode(zr, zi, Delta, pt);
  double expected = sqrt(K * N / 12.0)  / (double) Delta * norm_pt(pt);
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 100; ++r) {
    double Ar[K][N/2], Ai[K][N/2];
    double Azr[N/2], Azi[N/2];
    double pt_Az_raw[N], pt_Az[N];
    double e_Az[N];
    
    // get pt_Az_raw
    set_random_matrix(Ar, Ai);
    matrix_vector_product(zr, zi, Ar, Ai, Azr, Azi);
    encode(Azr, Azi, DeltaSq, pt_Az_raw);

    // get pt_Az
    matrix_vector_product(pt, Ar, Ai, pt_Az);

    // get e_Az
    sub_pt(pt_Az_raw, pt_Az, e_Az);

    // sanity check
    double er[N/2], ei[N/2];
    decode(e_Az, DeltaSq, er, ei);
    //print_pt("er (for sanity check)", er);

    double measured = norm(er, ei);
    std::cout << measured << std::endl;
    //double bound = K * expected;
    //std::cout << "Measured : " << measured << std::endl;
    //std::cout << "Expected : " << expected << std::endl;
    //std::cout << "Bound : " << bound << std::endl;
    //std::cout << "Measured / Expected : " << (measured / expected) << std::endl;
  } 
}
