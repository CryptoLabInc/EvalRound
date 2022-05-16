#pragma once

#include "HEAAN.h"
#include "DFT.h"

template<int N>
void rot( const int s[N], int s_rot[N] ){
	for(int i=0; i<N; i++){
		int q = (5*i)/N;
		int r = (5*i)%N;
		s_rot[r] = (q%2==0)?s[i]:-s[i];
	}
}

template<int N>
void rot( const int s[N], int s_rot[N], const int r){
	if(r == 0) {
		for(int i = 0; i < N; ++i)
			s_rot[i] = s[i];
			return;
	}

	int s_tmp[N];
	rot<N>(s, s_rot);
	for(int i = 1; i < r; ++i) {
		for(int j = 0; j < N; ++j) {
			s_tmp[j] = s_rot[j];
		}
		rot<N>(s_tmp, s_rot);
	}
}

template<int LOGQ, int N>
void rot( const R_Q<LOGQ,N>& pt, 
				R_Q<LOGQ,N>& pt_rot ){
	for(int i=0; i<N; i++){
		int q = (5*i)/N;
		int r = (5*i)%N;
		pt_rot[r] = pt[i];
		if(q%2==1) pt_rot[r].negate();
	}
}

template<int LOGQ, int N>
void rot_ct( const R_Q_square<  LOGQ,N>& ct,
			 const R_Q_square<2*LOGQ,N>& rkey,
			   	   R_Q_square<  LOGQ,N>&ct_rot ){
	rot<LOGQ,N>(ct[0],ct_rot[0]);
	rot<LOGQ,N>(ct[1],ct_rot[1]);
	HEAAN<LOGQ,N>::ks(rkey,ct_rot);
}

//-------------------------------------------------------------------------------
// rotate ct off times
//-------------------------------------------------------------------------------
template<int LOGQ, int N>
void rot_ct(const R_Q_square<  LOGQ, N>& ct, int off,
	        const R_Q_square<2*LOGQ, N>& rkey,
	              R_Q_square<  LOGQ, N>& ct_rot) {
	R_Q<LOGQ, N> temp; ct_rot = ct;
	for (int i = 0; i < off; i++) {
		rot<LOGQ, N>(ct_rot[0], temp); ct_rot[0] = temp;
		rot<LOGQ, N>(ct_rot[1], temp); ct_rot[1] = temp;
	}
	HEAAN<LOGQ, N>::ks(rkey, ct_rot);
}

template<int N>
void conj( const int s[N], int s_conj[N] ){
	s_conj[0] = s[0];
	for(int i=1; i<N; i++) s_conj[i] = -s[N-i];
}

template<int LOGQ, int N>
void conj( const R_Q<LOGQ,N>& pt,
				 R_Q<LOGQ,N>& pt_conj ){
	pt_conj[0] = pt[0];
	for(int i=1; i<N; i++){
		pt_conj[i] = pt[N-i];
		pt_conj[i].negate();
	}
}

template<int LOGQ, int N>
void conj( const R_Q_square<  LOGQ,N>& ct,
		   const R_Q_square<2*LOGQ,N>& ckey,
				 R_Q_square<  LOGQ,N>& ct_conj ){
	conj<LOGQ,N>(ct[0], ct_conj[0]);
	conj<LOGQ,N>(ct[1], ct_conj[1]);
	HEAAN<LOGQ,N>::ks(ckey,ct_conj);
}

template<int N>
void mat_to_vecs( const double A[N][N],
						double V[N][N]){
	for(int i=0; i<N; i++)
	for(int j=0; j<N; j++)
		V[i][j] = A[j][(i+j)%N];
}

template<int LOGQ, int LOGN, int LOGDELTA>
void linear_transform( const double Ar[1 << (LOGN - 1)][1 << (LOGN - 1)],
					   const double Ai[1 << (LOGN - 1)][1 << (LOGN - 1)],
					   const R_Q_square<  LOGQ, 1 << LOGN>& ct,
					   const R_Q_square<2*LOGQ, 1 << LOGN>& rkey,
							 R_Q_square<  LOGQ, 1 << LOGN>& Act ){
	const int N = 1 << LOGN;

	double vr[N/2][N/2]; mat_to_vecs<N/2>(Ar,vr);
	double vi[N/2][N/2]; mat_to_vecs<N/2>(Ai,vi);
	Act.setzero(); R_Q_square<LOGQ,N> ct_rot = ct;
	for(int i=0; i<N/2; i++){
		R_Q<LOGQ,N> pt;
		encode<LOGQ, LOGN>(vr[i],vi[i], 1ULL<<LOGDELTA, pt);

		

		R_Q_square<LOGQ,N> ct_rot_pt = ct_rot;
		ct_rot_pt *= pt;

		Act += ct_rot_pt; R_Q_square<LOGQ,N> temp; 
		rot_ct<LOGQ,N>(ct_rot, rkey, temp); // BUG
		ct_rot = temp; 
	}
}


//-----------------------------------------------------------------------------------
// A : sparse diagonal matrix
// Act: lineartransform(A,ct) such that
// decode(dec(Act,s,Q),Delta^2) = A. decode(dec(ct,s,Q),Delta)
//-----------------------------------------------------------------------------------
template< int LOGQ, int LOGN, int LOGDELTA, int S >
void linear_transform( const SparseDiagonal<1 << (LOGN - 1),S>& Ar,
					   const SparseDiagonal<1 << (LOGN - 1),S>& Ai,
					   const R_Q_square<  LOGQ, 1 << LOGN>& ct,
					   const R_Q_square<2*LOGQ, 1 << LOGN>* rkey[S],
							 R_Q_square<  LOGQ, 1 << LOGN>& Act ){
	const int N = 1 << LOGN;
	Act.setzero();
	for (int s = 0; s < S; s++) {
		if (Ar.zero[s] == false || Ai.zero[s] == false) {
			R_Q<LOGQ, N> pt;
			encode<LOGQ, LOGN>(Ar.vec[s],
							Ai.vec[s], 1ULL << LOGDELTA, pt);
			R_Q_square<LOGQ, N> ct_rot;
			if (Ar.off[s] != 0)
				rot_ct<LOGQ, N>(ct, Ar.off[s],*(rkey[s]), ct_rot);
			ct_rot *= pt;
			Act += ct_rot;
		}
	}
}


template< int LOGQ, int LOGN, int LOGDELTA, int S >
void linear_transform( const SparseDiagonal<1 << (LOGN - 1),S>& Ar,
					   const SparseDiagonal<1 << (LOGN - 1),S>& Ai,
					   const R_Q_square<  LOGQ, 1 << LOGN>& ct,
						const int skey[1 << LOGN],
							 R_Q_square<  LOGQ, 1 << LOGN>& Act ){
	const int N = 1 << LOGN;
	Act.setzero();
	for (int s = 0; s < S; s++) {
		if (Ar.zero[s] == false || Ai.zero[s] == false) {
			R_Q<LOGQ, N> pt;
			encode<LOGQ, LOGN>(Ar.vec[s],
							Ai.vec[s], 1ULL << LOGDELTA, pt);
			R_Q_square<LOGQ, N> ct_rot(ct);
			if (Ar.off[s] != 0) {
				R_Q_square<2*LOGQ, 1 << LOGN> rkey;
				int skey_rot[N];
				rot<N>(skey, skey_rot, Ar.off[s]);
				HEAAN<LOGQ,N>::swkgen(skey_rot ,skey, rkey);
				rot_ct<LOGQ, N>(ct, Ar.off[s], rkey, ct_rot);
			}
			ct_rot *= pt;
			Act += ct_rot;
		}
	}
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void serial_linear_transform( const SparseDiagonal<1 << (LOGN - 1),S> Ar[D],
					   	const SparseDiagonal<1 << (LOGN - 1),S> Ai[D],
					   	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
						const int skey[1 << LOGN],
						R_Q_square<  LOGQ, 1 << LOGN>& Act){
	R_Q_square<LOGQ, 1 << LOGN> ct_temp;
	Act = ct;
	for(int d = 0; d < D; ++d) {
		ct_temp = Act;
		linear_transform<LOGQ, LOGN, LOGDELTA, S>(Ar[d], Ai[d], ct_temp, skey, Act);
	} 
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped_serial_linear_transform( const SparseDiagonal<1 << (LOGN - 1),S> Ar[D],
					   	const SparseDiagonal<1 << (LOGN - 1),S> Ai[D],
					   	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
						const int skey[1 << LOGN],
						R_Q_square<  LOGQ, 1 << LOGN>& Act){
	// TODO(ksh) : implement this for all D
	assert(D % 2 == 0);
	SparseDiagonal<1 << (LOGN - 1), S*S> Br[D/2];
	SparseDiagonal<1 << (LOGN - 1), S*S> Bi[D/2];

	for(int i = 0; i < D/2; ++i) {
		MatMul(Ar[2*i], Ai[2*i], Ar[2*i+1], Ai[2*i+1], Br[i], Bi[i]);
	}
	serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D/2>(Br, Bi, ct, skey, Act);
}