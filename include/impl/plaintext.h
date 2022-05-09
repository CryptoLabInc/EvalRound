#pragma once

#include <iostream>
#include <string>

#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_matrix.h"

template <int LOGQ, int N>
void rotate(const R_Q<LOGQ,N>& pt, 
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

template <int LOGQ, int LOGN, int N, int K>
void matrix_vector_product(
    const R_Q<LOGQ, N> &pt,
    const Message<LOGN> A[k],
    R_Q<LOGQ, N> &pt_Az){
	R_Q<LOGQ, N> pt_v, pt_rot, pt_conv;
	pt_Az.setzero();
	for(int k = 0; k < K; ++k) {
		encode<LOGQ, LOGN>(A[k].r, A[k].i, Delta, pt_v);
		rotate<LOGQ, N>(pt, pt_rot, k);
		conv<LOGQ, N>(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}

template <int LOGQ, int LOGN, int N, int K>
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
		rotate(pt, pt_rot, Ar.off[k]);
		conv<LOGQ, N>(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}
