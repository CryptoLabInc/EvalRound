#pragma once

#include <cmath>

#include "matrix/DFT.h"
#include "arith/R_Q_square.h"

template< int LOGQ, int LOGN >
void encode( const double zr[1 << (LOGN - 1)], 
			 const double zi[1 << (LOGN - 1)], uint64_t Delta, R_Q<LOGQ, 1<< LOGN >& pt){
	const int N = 1 << LOGN;
	double m[N]; ifft<LOGN>(zr,zi,m);
	for(int i=0; i<N; i++){
		pt[i].setzero();
		pt[i][0]= (uint64_t) std::round(abs(m[i])*Delta);
		if(m[i]<0)
			pt[i].negate();
	}
}

template< int LOGQ, int LOGN >
void decode( const R_Q<LOGQ,1 << LOGN>& pt, uint64_t Delta, double zr[1 << (LOGN - 1)], 
													double zi[1 << (LOGN - 1)]){
	const int N = 1 << LOGN;
	double m[N];
	uint64_t alpha=1; alpha<<=32;
	double beta = (double)alpha; beta=beta*beta;
	for(int i=0; i<N; i++){
		Z_Q<LOGQ> abs_pti(pt[i]);
		if(pt[i].is_bigger_than_halfQ())
			abs_pti.negate();
		m[i]=0; double power=1;
		for(int j=0; j<(LOGQ+63)/64; j++, power*=beta)
			if(abs_pti[j]!=0)
				m[i]+=abs_pti[j]*power/Delta;
		if(pt[i].is_bigger_than_halfQ())
			m[i]=-m[i];
	}
	fft<LOGN>(m,zr,zi);
}

template< int LOGQ, int LOGN >
void decode_log( const R_Q<LOGQ, 1<<LOGN >& pt, int LOGDELTA, double zr[1 << (LOGN - 1)], 
													double zi[1 << (LOGN - 1)]){
	const int N = 1 << LOGN;
	double m[N];
	uint64_t alpha=1; alpha<<=32;
	double beta = (double)alpha; beta=beta*beta;
	for(int i=0; i<N; i++){
		Z_Q<LOGQ> abs_pti(pt[i]);
		if(pt[i].is_bigger_than_halfQ())
			abs_pti.negate();
		int d =abs_pti.max_valid_digit();
		int scale_diff = 64*d - LOGDELTA;
		m[i]=0;
		for(int j=0; j<= d; j++)
			m[i] = m[i] / beta + abs_pti[j];
		m[i] = ldexp(m[i], scale_diff);
		if(pt[i].is_bigger_than_halfQ())
			m[i]=-m[i];
	}
	fft<LOGN>(m,zr,zi);
}