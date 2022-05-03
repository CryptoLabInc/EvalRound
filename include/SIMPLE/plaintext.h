#pragma once

#include <iostream>
#include <string>

#include "HEAAN/DFT.h"
#include "define.h"

void encode( const double zr[N/2],
			 const double zi[N/2], uint64_t Delta, double pt[N]){
	double m[N]; ifft<LOGN>(zr,zi,m);
	double Delta_double = (double) Delta;
	for(int i=0; i<N; i++){
		pt[i] = (double) round(m[i]*Delta_double);
	}
}

void encode_raw( const double zr[N/2], 
			 const double zi[N/2], uint64_t Delta, double pt[N]){
	double m[N]; ifft<LOGN>(zr,zi,m);
	double Delta_double = (double) Delta;
	for(int i=0; i<N; i++){
		pt[i] = m[i]*Delta_double;
	}
}

void decode( const double pt[N], uint64_t Delta, double zr[N/2], 
													double zi[N/2]){
	double m[N];
	double Delta_double = (double) Delta;
    for(int i = 0; i < N; ++i) {
        m[i] = pt[i] / Delta_double;
    }
	fft<LOGN>(m,zr,zi);
}

void conv(const double pt1[N], const double pt2[N], double pt3[N]){
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int k=0;   k<=i; k++) sum += pt1[k]*pt2[i-k];
		for(int k=i+1; k< N; k++) sum -= pt1[k]*pt2[i+N-k];
		pt3[i]=sum;
	}
}

void set_zero(double pt[N]) {
	for(int i = 0; i <N; ++i) {
		pt[i] = 0;
	}
}

void add_pt(const double pt1[N], const double pt2[N], double pt3[N]){
	for(int i=0; i<N; i++){
		pt3[i] = pt1[i] + pt2[i];
	}
}

void sub_pt(const double pt1[N], const double pt2[N], double pt3[N]){
	for(int i=0; i<N; i++){
		pt3[i] = pt1[i] - pt2[i];
	}
}

double square_sum_pt(const double pt[N]) {
	double sum = 0;
	for(int i = 0; i < N; ++i) {
		sum += pt[i]* pt[i];
	}
	return sum;
}

double norm_pt(const double pt[N]){
    return sqrt(square_sum_pt(pt));
}

void rotate_once(const double pt[N], double pt_rot[N]){
	for(int i=0; i<N; i++){
		int q = (5*i)/N;
		int r = (5*i)%N;
		pt_rot[r] = pt[i];
		if(q%2==1) pt_rot[r] = -pt_rot[r];
	}
}

void rotate_pt(const double pt[N], double pt_rot[N], int r){
	double pt_tmp[N];
	for(int i = 0; i < N; ++i) {
		pt_rot[i] = pt[i];
	}
	for(int rot = 0; rot < r; ++rot) {
		for(int i = 0; i < N; ++i) {
			pt_tmp[i] = pt_rot[i];
		}
		rotate_once(pt_tmp, pt_rot);
	}
}

void matrix_vector_product(
    const double pt[N],
    const double Ar[K][N/2], const double Ai[K][N/2],
    double pt_Az[N]){
	double pt_v[N], pt_rot[N], pt_conv[N];
	set_zero(pt_Az);
	for(int k = 0; k < K; ++k) {
		encode(Ar[k], Ai[k], Delta, pt_v);
		rotate_pt(pt, pt_rot, k);
		conv(pt_v, pt_rot, pt_conv);
		add_pt(pt_Az, pt_conv, pt_Az);
	}
}

void matrix_vector_product(
    const double pt[N],
    SparseDiagonal<(1<<(LOGN-1)),3> Ar,
	SparseDiagonal<(1<<(LOGN-1)),3> Ai,
    double pt_Az[N]) {
	double pt_v[N], pt_rot[N], pt_conv[N];
	set_zero(pt_Az);
	for(int k = 0; k < 3; ++k) {
		assert(Ar.off[k] == Ai.off[k]);
		encode(Ar.vec[k], Ai.vec[k], Delta, pt_v);
		rotate_pt(pt, pt_rot, Ar.off[k]);
		conv(pt_v, pt_rot, pt_conv);
		add_pt(pt_Az, pt_conv, pt_Az);
	}
}

void print_pt(const std::string name, const double pt[N]){
	std::cout << "Plaintext " << name << std::endl;
	for(int i=0; i <std::min(N, 10); ++i) {
		std::cout << pt[i] << " ";
	}
	std::cout << std::endl;
}