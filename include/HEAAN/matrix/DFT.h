#pragma once

#include "matrix.h"
#include "util/util.h"

#include <math.h>
#include <stdio.h>
#include <assert.h>

template<int N>
void dft( const double m[N], double zr[N/2],
							 double zi[N/2]){
	int fivei=1;
	for(int i=0;i<N/2;i++,fivei=(fivei*5)%(2*N)){
		int fivei_k=0; zr[i]=0; zi[i]=0;
		for(int k=0;k<N;k++,fivei_k=(fivei_k+fivei)%(2*N)){
			zr[i]+=m[k]*cos(PI/N*fivei_k);
			zi[i]+=m[k]*sin(PI/N*fivei_k);
		}
	}
}

template<int N>
void idft( const double zr[N/2],
		   const double zi[N/2], double m[N]){
	for(int k=0;k<N/2;k++){
		double sum1 = 0, sum2 = 0;
		int fivei_k=k;
		for(int i=0;i<N/2;i++,fivei_k=(fivei_k*5)%(2*N)){
			sum1+=cos(PI/N*fivei_k)*zr[i]+sin(PI/N*fivei_k)*zi[i];
			sum2+=cos(PI/N*fivei_k)*zi[i]-sin(PI/N*fivei_k)*zr[i];
		}
		m[k    ] =2./N * sum1;
		m[k+N/2] =2./N * sum2;
	}
}

//------------------------------------------------------------------------------------------
//
// When N=16, U0NR matrix is decomposed into the product of three 8x8 complex matrices.
// U0NR = A[0] x A[1] x A[2]
// Each A[i] is a sparse diagonal matrix with three diagonal vectors.
//------------------------------------------------------------------------------------------
template< int LOGN >
void splitU0NR( SparseDiagonal<(1<<(LOGN-1)),3> Ar[LOGN-1],
	            SparseDiagonal<(1<<(LOGN-1)),3> Ai[LOGN-1]){
	for (int n = 0; n < LOGN - 1; n++) {
		Ar[n].setzero(); int N = 1 << LOGN;
		Ai[n].setzero(); int M = 1 <<(LOGN-n-2);
		Ar[n].off[0] =    0; Ai[n].off[0] =    0;
		Ar[n].off[1] =    M; Ai[n].off[1] =    M;
		Ar[n].off[2] =N/2-M; Ai[n].off[2] =N/2-M;

		for (int k = 0; k < (1 << n); k++)
			for (int i = 0, fivei = (1 << n); i < M; i++, fivei = (fivei * 5) % (2 * N)) {
				Ar[n].set_element(2 * M*k +     i, 2 * M*k +     i, 1);
				Ar[n].set_element(2 * M*k + M + i, 2 * M*k +     i, 1);
				Ar[n].set_element(2 * M*k     + i, 2 * M*k + M + i, cos(PI/N*fivei));
				Ai[n].set_element(2 * M*k     + i, 2 * M*k + M + i, sin(PI/N*fivei));
				Ar[n].set_element(2 * M*k + M + i, 2 * M*k + M + i,-cos(PI/N*fivei));
				Ai[n].set_element(2 * M*k + M + i, 2 * M*k + M + i,-sin(PI/N*fivei));
			}
	}
}

template<int LOGN>
void fftNR(
  const double m[1 << LOGN],
  double zr[1 << (LOGN - 1)],
  double zi[1 << (LOGN - 1)]){
	const int N = 1 << LOGN;
	static bool is_init = false;
	static SparseDiagonal<N/2,3> U0r[LOGN-1];
    static SparseDiagonal<N/2,3> U0i[LOGN-1];
	if(!is_init) {
  		splitU0NR<LOGN>(U0r, U0i);
		is_init = true;
	}

	for(int i = 0; i < N/2; ++i) {
		zr[i] = m[i];
		zi[i] = m[i + N/2];
	}

	double Umr[N/2], Umi[N/2];
	for(int d = LOGN - 2; d >= 0; --d) {
		for(int i = 0; i < N/2; ++i) {
			Umr[i] = zr[i];
			Umi[i] = zi[i];
		}
		matrix_vector_product<LOGN>(Umr, Umi, U0r[d], U0i[d], zr, zi);
	}
}


template<int LOGN>
void ifftNR( const double zr[1 << (LOGN - 1)],
		   const double zi[1 << (LOGN - 1)], double m[1 << LOGN]){
	const int N = 1 << LOGN;
	static bool is_init = false;
	static SparseDiagonal<N/2,3> U0r[LOGN-1];
    static SparseDiagonal<N/2,3> U0i[LOGN-1];

	if(!is_init) {
		splitU0NR<LOGN>(U0r, U0i);
    	for (int n = 0; n < LOGN - 1; n++) {
			U0r[n].transpose();
			U0i[n].transpose();
			U0i[n].negate();
		}
		is_init = true;
	}

	for(int i = 0; i < N/2; ++i) {
		m[i] = zr[i];
		m[i + N/2] = zi[i];
	}

	double Uzr[N/2], Uzi[N/2];
	for(int d = 0; d < LOGN - 1; ++d) {
		for(int i = 0; i < N/2; ++i) {
			Uzr[i] = m[i];
			Uzi[i] = m[i + N/2];
		}
		matrix_vector_product<LOGN>(Uzr, Uzi, U0r[d], U0i[d], m, m + N/2);
	}

	for(int i = 0; i < N; ++i) {
		m[i] = m[i] * 2 / (double) N;
	}
}

template<int LOGN>
void fft(
  const double m[1 << LOGN],
  double zr[1 << (LOGN - 1)],
  double zi[1 << (LOGN - 1)]){
	const int N = 1 << LOGN;
	double mNR[N];
	for(uint32_t i = 0; i < N/2; ++i) {
		uint32_t i_NR = bitReverse(i, LOGN - 1);
		mNR[i_NR] = m[i];
		mNR[i_NR+N/2] = m[i+N/2];
	}
	fftNR<LOGN>(mNR, zr, zi);
}

template<int LOGN>
void ifft( const double zr[1 << (LOGN - 1)],
		   const double zi[1 << (LOGN - 1)], double m[1 << LOGN]){
	const int N = 1 << LOGN;
	double mNR[N];
	ifftNR<LOGN>(zr, zi, mNR);
	for(uint32_t i = 0; i < N/2; ++i) {
		uint32_t i_NR = bitReverse(i, LOGN - 1);
		m[i] = mNR[i_NR];
		m[i+N/2] = mNR[i_NR+N/2];
	}
}