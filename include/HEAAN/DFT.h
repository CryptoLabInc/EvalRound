#pragma once
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "matrix.h"

#define PI 3.1415926535897932384626433

#include <iostream>

template<int N>
void dft_orig( const double m[N], double zr[N/2],
							 double zi[N/2]){
	for(int i=0;i<N/2;i++){
		double sumr = 0, sumi = 0;
		int idxi_k=0;
		for(int k=0;k<N;k++,idxi_k=(idxi_k+2*i+1)%(2*N)){
			sumr+=m[k]*cos(PI/N*idxi_k);
			sumi+=m[k]*sin(PI/N*idxi_k);
		}
		zr[i] = sumr;
		zi[i] = sumi;
	}
}

template<int N>
void fft_orig( const double m[N], double zr[N/2],
							 double zi[N/2]){
	if(N == 2){
		zr[0] = m[0];
		zi[0] = m[1];
		return;
	}
	double e[N/2], o[N/2];
	double e_zr[N/4], e_zi[N/4], o_zr[N/4], o_zi[N/4];
	for(int i=0; i < N/2; ++i) {
		e[i] = m[2*i];
		o[i] = m[2*i + 1];
	}
	fft_orig<N/2>(e, e_zr, e_zi);
	fft_orig<N/2>(o, o_zr, o_zi);
	for(int i = 0; i < N/4; ++i){
		const double rot = PI/N*(2*i+1);
		zr[i] = e_zr[i] + o_zr[i] * cos(rot) - o_zi[i] * sin(rot);
		zi[i] = e_zi[i] + o_zr[i] * sin(rot) + o_zi[i] * cos(rot);
	}
	for(int i = N/4; i < N/2; ++i){
		const double rot = PI/N*(2*i+1);
		int j = N/2 - 1 - i; // conjugate index
		zr[i] = e_zr[j] + o_zr[j] * cos(rot) + o_zi[j] * sin(rot);
		zi[i] = -e_zi[j] + o_zr[j] * sin(rot) - o_zi[j] * cos(rot);
	}
}


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
	double zr_orig[N/2], zi_orig[N/2];
	fft_orig<N>(m, zr_orig, zi_orig);
	int fivei = 1;
	for(int i=0;i<N/2;i++,fivei=(fivei*5)%(2*N)){
		if(fivei < N) {
			zr[i] = zr_orig[fivei / 2];
			zi[i] = zi_orig[fivei / 2];
		} else {
			zr[i] = zr_orig[(2*N - fivei) / 2];
			zi[i] = -zi_orig[(2*N - fivei) / 2];
		}
	}
}

template<int N>
void idft_orig( const double zr[N/2],
		   const double zi[N/2], double m[N]){
	for(int k=0;k<N/2;k++){
		double sum_real = 0, sum_imag = 0;
		m[k]=m[k+N/2]=0; int idxi_k=k;
		for(int i=0;i<N/2;i++,idxi_k=(idxi_k+2*k)%(2*N)){
			sum_real +=cos(PI/N*idxi_k)*zr[i]+sin(PI/N*idxi_k)*zi[i];
			sum_imag +=cos(PI/N*idxi_k)*zi[i]-sin(PI/N*idxi_k)*zr[i];
		}
		m[k    ] = 2./N * sum_real;
		m[k+N/2] = 2./N * sum_imag;
	}
}

template<int N>
void ifft_orig(const double zr[N/2],
		const double zi[N/2], double m[N]){
	if(N == 2){
		m[0] = zr[0];
		m[1] = zi[0];
		return;
	}
	double e_zr[N/4], e_zi[N/4], o_zr[N/4], o_zi[N/4];
	double e[N/2], o[N/2];
	for(int i=0; i < N/4; ++i) {
		e_zr[i] = zr[2*i];
		e_zi[i] = zi[2*i];
		o_zr[i] = zr[2*i+1];
		o_zi[i] = zi[2*i+1];
	}
	ifft_orig<N/2>(e_zr, e_zi, e);
	ifft_orig<N/2>(o_zr, o_zi, o);
	for(int k = 0; k < N/4; ++k){
		const double rot_even = -PI/N*k;
		const double rot_odd = PI/N*k;
		m[k] = e[k] * cos(rot_even) + e[k+N/4] * sin(rot_even);
		m[k] += o[k] * cos(rot_odd) + o[k+N/4] * sin(rot_odd);
		m[k+N/2] = - e[k] * sin(rot_even) + e[k+N/4] * cos(rot_even);
		m[k+N/2] += - o[k] * sin(rot_odd) + o[k+N/4] * cos(rot_odd);
		m[k    ]*=0.5;
		m[k+N/2]*=0.5;
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

