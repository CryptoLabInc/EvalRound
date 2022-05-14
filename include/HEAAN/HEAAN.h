#pragma once
#include <math.h>
#include <stdlib.h>
#include "R_Q_square.h"
#include "DFT.h"

template< int LOGQ, int N>
struct HEAAN{
	static void keygen( int h, int s[N] );
	static void enc( const R_Q       <LOGQ,N>& pt, const int s[N], R_Q_square<LOGQ,N>& ct );
	static void dec( const R_Q_square<LOGQ,N>& ct, const int s[N], R_Q       <LOGQ,N>& pt );
	static void swkgen( const int s_fr[N], const int s_to[N], R_Q_square<2*LOGQ,N>& swk );
	static void ks    ( const R_Q_square<2*LOGQ,N>& swk, R_Q_square<LOGQ,N>& ct );
};

//-------------------------------------------------------------------------------------
// implementations
//-------------------------------------------------------------------------------------
template< int LOGQ, int N >
void HEAAN<LOGQ,N>::keygen( int h, int s[N] ){
	for(int i=0;i<N;i++) s[i]=0;
	for(int i=0;i<h;i++){
		int j=rand()%N;
		while(s[j]!=0)
			j=rand()%N;
		s[j]=(rand()%2==0)?1:-1;
	}
}


template< int LOGQ, int N >
void HEAAN<LOGQ,N>::swkgen( const int s_fr[N],
						    const int s_to[N],
						    R_Q_square<2*LOGQ,N>& swk){
	Z_Q<  LOGQ> one; one.setzero(); one[0]=1;
	Z_Q<2*LOGQ> P;
	shift_left<LOGQ,LOGQ*2>(one,P);
	
	R_Q<2*LOGQ,N> Ps_fr;
	for(int i=0; i<N; i++)
		if  (s_fr[i]== 0)  Ps_fr[i].setzero();
		else{
			int abs=(s_fr[i]>0)? s_fr[i]:-s_fr[i];
			Ps_fr[i]=P; abs--;
			while(abs>0){
				Ps_fr[i]+=P; abs--;
			}
			if(s_fr[i]<0) Ps_fr[i].negate();
		}
	HEAAN<2*LOGQ,N>::enc(Ps_fr,s_to,swk);
}


uint64_t rand_uint64(void){
	uint64_t r =0;
	for( int i=0; i<64; i+=15 /*30*/){
		  r = r*((uint64_t)RAND_MAX + 1) + rand();
  }
  return r;
}

template< int LOGQ, int N >
void HEAAN<LOGQ,N>::enc(const R_Q<LOGQ,N>& pt,
						const int s[N],
							  R_Q_square<LOGQ,N>& ct){
	for(int i=0; i<N; i++){
		for(int j=0; j<(LOGQ+63)/64; j++)
			ct[1][i][j] = rand_uint64();
		ct[1][i].remove_clutter();
	}
	R_Q<LOGQ,N> e;
	for(int i=0; i<N; i++){
		e[i].setzero();
		e[i][0] = rand()%9;
		if(rand()%2==1)
			e[i].negate();
	}
	e    +=pt;
	ct[0] =e;
	conv<LOGQ,N> (s,ct[1],e);
	ct[0]-=e;
}



template< int LOGQ, int N >
void HEAAN<LOGQ,N>::dec(const R_Q_square <LOGQ,N>& ct,
						const int s[N],
							  R_Q<LOGQ,N>& pt){
	conv<LOGQ,N> (s,ct[1],pt);
	pt+=ct[0];
}


template< int LOGQfr, int LOGQto, int N>
void RS( const R_Q_square<LOGQfr,N>& ct_fr,
			   R_Q_square<LOGQto,N>& ct_to){
	for(int i=0; i<N; i++){
		shift_right<LOGQfr, LOGQto>(ct_fr[0][i],ct_to[0][i]);
		shift_right<LOGQfr, LOGQto>(ct_fr[1][i],ct_to[1][i]);
	}
}

template< int LOGQfr, int LOGQto, int N >
void mod_raise(const R_Q<LOGQfr,N>& ptfr,
					 R_Q<LOGQto,N>& ptto){
	for(int i=0;i<N             ;i++)
	for(int j=0;j<(LOGQto+63)/64;j++)
		if(j<(LOGQfr+63)/64) ptto[i][j]=ptfr[i][j];
		else                 ptto[i][j]=0;
}

template< int LOGQfr, int LOGQto, int N >
void mod_raise(const R_Q_square<LOGQfr,N>& ctfr,
					 R_Q_square<LOGQto,N>& ctto){
	R_Q<LOGQto, N> pt0, pt1;
    mod_raise<LOGQfr, LOGQto, N>(ctfr[0], pt0);
	ctto[0] = pt0;
	mod_raise<LOGQfr, LOGQto, N>(ctfr[1], pt1);
	ctto[1] = pt1;
}

template< int LOGQ, int N >
void HEAAN<LOGQ,N>::ks(const R_Q_square <2*LOGQ,N>& swk,
							 R_Q_square <  LOGQ,N>& ct){
	R_Q_square<2*LOGQ,N> ct1_swk;
	ct1_swk = swk;
	R_Q<2*LOGQ,N> ct1_;
	for(int i=0;i<N             ;i++)
	for(int j=0;j<(2*LOGQ+63)/64;j++)
		if(j<(LOGQ+63)/64) ct1_[i][j]=ct[1][i][j];
		else               ct1_[i][j]=0;
	ct1_swk*=ct1_;
	R_Q_square<LOGQ,N> ct1_swk_rs;
	RS<LOGQ*2,LOGQ,N>(ct1_swk,ct1_swk_rs);
	for(int i=0; i<N; i++) ct[1][i].setzero();
	ct+=ct1_swk_rs;
}

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
	int scale_diff = 64 * ((LOGQ+63)/64 - 1) - LOGDELTA;
	for(int i=0; i<N; i++){
		Z_Q<LOGQ> abs_pti(pt[i]);
		if(pt[i].is_bigger_than_halfQ())
			abs_pti.negate();
		m[i]=0;
		for(int j=0; j<(LOGQ+63)/64; j++)
			m[i] = m[i] / beta + abs_pti[j];
		m[i] = ldexp(m[i], scale_diff);
		if(pt[i].is_bigger_than_halfQ())
			m[i]=-m[i];
	}
	fft<LOGN>(m,zr,zi);
}

template<int N>
void conv( const int s1[N],
		   const int s2[N], int s3[N]){
	for(int i=0; i<N; i++){
		s3[i]=0;
		for(int k=0;   k<=i; k++) s3[i] += s1[k]*s2[i-k];
		for(int k=i+1; k< N; k++) s3[i] -= s1[k]*s2[i+N-k];
	}
}
template<int LOGQ, int N>
void Mul(const R_Q_square<  LOGQ,N>& ct1,
		 const R_Q_square<  LOGQ,N>& ct2,
		 const R_Q_square<2*LOGQ,N>& evk,
			   R_Q_square<  LOGQ,N>& ct3){
	R_Q<LOGQ,N> d2; conv<LOGQ,N>(ct1[1],ct2[1],d2);
	R_Q<2*LOGQ,N> d2_;
	mod_raise<LOGQ, 2*LOGQ,N> (d2,d2_);
	R_Q_square<2*LOGQ,N> evk_d2(evk);
	evk_d2*=d2_;
	RS<2*LOGQ, LOGQ, N> (evk_d2, ct3);
	R_Q_square<LOGQ,N> d0d1;
	conv<LOGQ,N>(ct1[0],ct2[1],d0d1[0]);
	conv<LOGQ,N>(ct1[1],ct2[0],d0d1[1]);
	d0d1[1]+=d0d1[0];
	conv<LOGQ,N>(ct1[0],ct2[0],d0d1[0]);
	ct3+=d0d1;
}

