#pragma once

#include "arith/R_Q_square.h"
#include "util/util.h"

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

template< int LOGQfr, int LOGQto, int N>
void RS( const R_Q_square<LOGQfr,N>& ct_fr,
			   R_Q_square<LOGQto,N>& ct_to){
	for(int i=0; i<N; i++){
		shift_right<LOGQfr, LOGQto>(ct_fr[0][i],ct_to[0][i]);
		shift_right<LOGQfr, LOGQto>(ct_fr[1][i],ct_to[1][i]);
	}
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

template<int LOGQ, int N>
void Mul(const R_Q_square<  LOGQ,N>& ct1,
		 const R_Q_square<  LOGQ,N>& ct2,
		 const R_Q_square<2*LOGQ,N>& evk,
			   R_Q_square<  LOGQ,N>& ct3){
	R_Q<LOGQ,N> d2; conv<LOGQ,N>(ct1[1],ct2[1],d2);
	R_Q<2*LOGQ,N> d2_;
	resize<LOGQ, 2*LOGQ,N> (d2,d2_);
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

