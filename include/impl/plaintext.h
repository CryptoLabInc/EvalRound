#pragma once

#include "HEAAN/linear_transform.h"

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

template <int LOGQ, int LOGN, int K>
void matrix_vector_product(
    const R_Q<LOGQ, 1 << LOGN> &pt,
    const Message<LOGN> A[K], uint64_t Delta,
    R_Q<LOGQ, 1 << LOGN> &pt_Az){
	R_Q<LOGQ, 1 << LOGN> pt_v, pt_rot, pt_conv;
	pt_Az.setzero();
	for(int k = 0; k < K; ++k) {
		encode<LOGQ, LOGN>(A[k].r, A[k].i, Delta, pt_v);
		rot<LOGQ, 1 << LOGN>(pt, pt_rot, k);
		conv(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}

template <int LOGQ, int LOGN>
void matrix_vector_product(
    const R_Q<LOGQ, 1 << LOGN> &pt,
    SparseDiagonal<1<<(LOGN-1),3> &Ar,
	SparseDiagonal<1<<(LOGN-1),3> &Ai, uint64_t Delta,
    R_Q<LOGQ, 1 << LOGN> &pt_Az) {
	R_Q<LOGQ, 1 << LOGN> pt_v, pt_rot, pt_conv;
	pt_Az.setzero();
	for(int k = 0; k < 3; ++k) {
		assert(Ar.off[k] == Ai.off[k]);
		encode<LOGQ, LOGN>(Ar.vec[k], Ai.vec[k], Delta, pt_v);
		rot(pt, pt_rot, Ar.off[k]);
		conv(pt_v, pt_rot, pt_conv);
		pt_Az += pt_conv;
	}
}
