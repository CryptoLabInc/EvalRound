#pragma once

#include <iostream>
#include <string>

#include "define.h"
#include "HEAAN/HEAAN.h"

double square_sum_pt(R_Q<LOGQ, N> pt) {
	// assume each Z_Qs < 2^64
	double sum = 0;
	for(int i = 0; i < N; ++i) {
		uint64_t val = pt.coeff[i][0];
		if(val >> 63)
			val = -val;
		sum += (double) val * (double) val;
	}
	return sum;
}

double norm_pt(const R_Q<LOGQ, N> pt){
    return sqrt(square_sum_pt(pt));
}

/*void rotate_once(const double pt[N], double pt_rot[N]){
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
    double pt_Az[N/2]){
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

void print_u(const std::string name, const double pt[N]){
	std::cout << "Plaintext " << name << std::endl;
	for(int i=0; i <std::min(N, 10); ++i) {
		std::cout << pt[i] << " ";
	}
	std::cout << std::endl;
}*/