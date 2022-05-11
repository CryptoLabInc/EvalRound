#pragma once

#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_matrix.h"
#include "HEAAN/advanced/conv.h"

#include <iostream>
#include <string>
#include <cmath>

template <int LOGQ, int LOGN>
void encode(const Message<LOGN> &z, uint64_t Delta, R_Q<LOGQ, 1<< LOGN >& pt){
	encode<LOGQ, LOGN>(z.r, z.i, Delta, pt);
}

template <int LOGQ, int LOGN>
void decode(const R_Q<LOGQ, 1 << LOGN>& pt, uint64_t Delta, Message<LOGN> &z){
	decode<LOGQ, LOGN>(pt, Delta, z.r, z.i);
}

template< int LOGQ, int LOGN >
void decode_log( const R_Q<LOGQ, 1<<LOGN >& pt, int LOGDELTA, Message<LOGN> &z){
	decode_log<LOGQ, LOGN>(pt, LOGDELTA, z.r, z.i);
}


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

template <int L, int LOGQ, int LOGN, int K>
void matrix_vector_product(
	const CONV<L, LOGQ, 1 << LOGN> &conv,
    const R_Q<LOGQ, 1 << LOGN> &pt,
    const Message<LOGN> A[K], uint64_t Delta,
    R_Q<LOGQ, 1 << LOGN> &pt_Az){
	R_Q<LOGQ, 1 << LOGN> pt_v, pt_rot, pt_conv;
	pt_Az.setzero();
	for(int k = 0; k < K; ++k) {
		encode<LOGQ, LOGN>(A[k].r, A[k].i, Delta, pt_v);
		rotate<LOGQ, 1 << LOGN>(pt, pt_rot, k);
		conv.conv(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}

template <int L, int LOGQ, int LOGN>
void matrix_vector_product(
	const CONV<L, LOGQ, 1 << LOGN> &conv,
    const R_Q<LOGQ, 1 << LOGN> &pt,
    SparseDiagonal<1<<(LOGN-1),3> &Ar,
	SparseDiagonal<1<<(LOGN-1),3> &Ai, uint64_t Delta,
    R_Q<LOGQ, 1 << LOGN> &pt_Az) {
	R_Q<LOGQ, 1 << LOGN> pt_v, pt_rot, pt_conv;
	pt_Az.setzero();
	for(int k = 0; k < 3; ++k) {
		assert(Ar.off[k] == Ai.off[k]);
		encode<LOGQ, LOGN>(Ar.vec[k], Ai.vec[k], Delta, pt_v);
		rotate(pt, pt_rot, Ar.off[k]);
		conv.conv(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}
