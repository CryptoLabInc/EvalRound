#pragma once

#include <iostream>
#include <string>

#include "define.h"
#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_matrix.h"

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

void rotate_pt(const R_Q<LOGQ,N>& pt, 
				R_Q<LOGQ,N>& pt_rot, const int r){
	if(r == 0) {
		pt_rot = pt;
		return;
	}
	R_Q<LOGQ, N> pt_tmp;
	rot<LOGQ, N>(pt, pt_rot);
	for(int i = 1; i < r; ++i) {
		pt_tmp = pt_rot;
		rot<LOGQ, N>(pt_tmp, pt_rot);
	}
}

void matrix_vector_product(
    const R_Q<LOGQ, N> &pt,
    const double Ar[K][N/2], const double Ai[K][N/2],
    R_Q<LOGQ, N> &pt_Az){
	R_Q<LOGQ, N> pt_v, pt_rot, pt_conv;
	pt_Az.setzero();
	for(int k = 0; k < K; ++k) {
		encode<LOGQ, LOGN>(Ar[k], Ai[k], Delta, pt_v);
		rotate_pt(pt, pt_rot, k);
		conv<LOGQ, N>(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}

void matrix_vector_product(
    const R_Q<LOGQ, N> &pt,
    SparseDiagonal<(1<<(LOGN-1)),3> Ar,
	SparseDiagonal<(1<<(LOGN-1)),3> Ai,
    R_Q<LOGQ, N> &pt_Az) {
	R_Q<LOGQ, N> pt_v, pt_rot, pt_conv;
	pt_Az.setzero();
	for(int k = 0; k < 3; ++k) {
		assert(Ar.off[k] == Ai.off[k]);
		encode<LOGQ, LOGN>(Ar.vec[k], Ai.vec[k], Delta, pt_v);
		rotate_pt(pt, pt_rot, Ar.off[k]);
		conv<LOGQ, N>(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}

/*void print_u(const std::string name, const double pt[N]){
	std::cout << "Plaintext " << name << std::endl;
	for(int i=0; i <std::min(N, 10); ++i) {
		std::cout << pt[i] << " ";
	}
	std::cout << std::endl;
}*/