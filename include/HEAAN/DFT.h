#pragma once
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "matrix.h"

#define PI 3.1415926535897932384626433

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
void fft( const double m[N], double zr[N/2],
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
		m[k]=m[k+N/2]=0; int fivei_k=k;
		for(int i=0;i<N/2;i++,fivei_k=(fivei_k*5)%(2*N)){
			m[k    ]+=cos(PI/N*fivei_k)*zr[i]+sin(PI/N*fivei_k)*zi[i];
			m[k+N/2]+=cos(PI/N*fivei_k)*zi[i]-sin(PI/N*fivei_k)*zr[i];
		}
		m[k    ]*=2./N;
		m[k+N/2]*=2./N;
	}
}

template<int N>
void ifft( const double zr[N/2],
		   const double zi[N/2], double m[N]){
	for(int k=0;k<N/2;k++){
		m[k]=m[k+N/2]=0; int fivei_k=k;
		for(int i=0;i<N/2;i++,fivei_k=(fivei_k*5)%(2*N)){
			m[k    ]+=cos(PI/N*fivei_k)*zr[i]+sin(PI/N*fivei_k)*zi[i];
			m[k+N/2]+=cos(PI/N*fivei_k)*zi[i]-sin(PI/N*fivei_k)*zr[i];
		}
		m[k    ]*=2./N;
		m[k+N/2]*=2./N;
	}
}

template<int N>
void get_U0( double U0r[N/2][N/2],
			 double U0i[N/2][N/2]){
	int fivei=1;
	for(int i=0;i<N/2;i++,fivei=(fivei*5)%(2*N)){
		int fivei_k=0;
		for(int k=0;k<N/2;k++,fivei_k=(fivei_k+fivei)%(2*N)){
			U0r[i][k]=cos(PI/N*fivei_k);
			U0i[i][k]=sin(PI/N*fivei_k);
		}
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

