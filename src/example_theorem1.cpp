
#include "BIG/test.h"

#include <iostream>

void norm_relationship_test() {
  double zr[N/2], zi[N/2];

  R_Q<LOGQ, N> pt;

  set_test_message(zr, zi);
  encode<LOGQ, N>(zr, zi, Delta, pt);

  double measured = norm_pt(pt) / norm(zr, zi);
  double expected = (double) Delta / sqrt(N/2);

  std::cout << "Measured : " << measured << std::endl;
  std::cout << "Expected : " << expected << std::endl;
}

void error_test(){
  double zr[N/2], zi[N/2];
  double zr_tilde[N/2], zi_tilde[N/2];
  double er[N/2], ei[N/2];
  R_Q<LOGQ, N> pt;

  // get pt
  set_random_message(zr, zi);
  encode<LOGQ, N>(zr, zi, Delta, pt);
  decode<LOGQ, N>(pt, Delta, zr_tilde, zi_tilde);
  sub(zr, zr_tilde, er);
  sub(zi, zi_tilde, ei);

  double measured = norm(er, ei);
  double expected = sqrt(N/2) / (double) Delta * sqrt(N/12.0);

  std::cout << "Measured : " << measured << std::endl;
  std::cout << "Expected : " << expected << std::endl;
}

void example(){
  // ||Az - \widetilde{Az} || ~= sqrt(KN / 12) * (1 / Delta_A) * ||z||
  double zr[N/2], zi[N/2];
  R_Q<LOGQ, N> pt;

  // get pt
  set_test_rounded_message(zr, zi);

  double expected = sqrt(K * N / 12.0) * norm(zr, zi) / (double) Delta;
  std::cout << expected << std::endl;
  std::cout << std::endl;

  for(int r = 0; r < 100; ++r) {
    double Ar[K][N/2], Ai[K][N/2];
    double Azr[N/2], Azi[N/2];
    R_Q<LOGQ, N> pt_Az;

    // somehow pt is changed ... reset pt
    encode<LOGQ, N>(zr, zi, Delta, pt);

    // get Az
    set_random_matrix(Ar, Ai);
    matrix_vector_product(zr, zi, Ar, Ai, Azr, Azi);

    // get pt_Az
    matrix_vector_product(pt, Ar, Ai, pt_Az);

    // decode
    double Azr_tilde[N/2], Azi_tilde[N/2];
    const int LOGQto = LOGQ - LOGDELTA;
    R_Q<LOGQto, N> pt_RS;
    for(int i = 0; i < N; ++i) {
		  shift_right<LOGQ, LOGQto>(pt_Az[i], pt_RS[i]);
	  }
    decode(pt_RS, Delta, Azr_tilde, Azi_tilde);
    //decode(pt_Az, Delta*Delta, Azr_tilde, Azi_tilde);
    /*R_Q<LOGQ, N> pt_Az_raw;
    encode(Azr, Azi, Delta*Delta, pt_Az_raw);
    decode(pt_Az_raw, Delta*Delta, Azr_tilde, Azi_tilde);*/

    //print("Azr", Azr);
    //print("Azr_tilde", Azr_tilde);

    // sub
    double er[N/2], ei[N/2];
    sub(Azr, Azr_tilde, er);
    sub(Azi, Azi_tilde, ei);

    double measured = norm(er, ei);
    std::cout << measured << " " << (measured / expected) << std::endl;

  }
}

int main()
{
  //error_test();
  example();
}
