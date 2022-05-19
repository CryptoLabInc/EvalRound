#pragma once

#include <math.h>
#include "HEAAN_matrix.h"
#include "HEAAN_poly.h"


template<int LOGQ, int LOGN>
void CoeffToSlot( const R_Q_square<  LOGQ,1<<LOGN>& ct,
				  const R_Q_square<2*LOGQ,1<<LOGN>& rkey,
				  const R_Q_square<2*LOGQ,1<<LOGN>& ckey,
						R_Q_square<  LOGQ,1<<LOGN> ct_[2]){
	const int N = 1 << LOGN;
	double U0r[N/2][N/2];
	double U0i[N/2][N/2]; get_U0<N>(U0r, U0i);
	double Ar[N/2][N/2];
	double Ai[N/2][N/2];

	for(int i=0; i<N/2; i++)
	for(int j=0; j<N/2; j++){
		Ar[i][j] =  U0r[j][i]/N;
		Ai[i][j] = -U0i[j][i]/N;
	}
	R_Q_square<LOGQ,N> ct1; linear_transform<LOGQ,LOGN,50>(Ar,Ai,ct,rkey,ct1);
	R_Q_square<LOGQ,N> ct2; conj(ct1,ckey,ct2);

	R_Q<LOGQ,N> pti; pti.setzero();
	pti[N/2][0] = 1; pti[N/2].negate();
	ct_[1] = ct1; ct_[1] -=ct2; ct_[1]*=pti;
	ct_[0] = ct1; ct_[0] +=ct2;
}

template<int LOGQ, int LOGN, int LOGDELTA, int G = 1>
void CoeffToSlot (	const R_Q_square<  LOGQ,1<<LOGN>& ct,
					const int s[1 << LOGN],
					R_Q_square<  LOGQ,1<<LOGN> ct_[2]){
	const int N = 1 << LOGN;
	static bool is_init = false;
	static SparseDiagonal<N/2,3> U0r[LOGN-1];
    static SparseDiagonal<N/2,3> U0i[LOGN-1];

	if(!is_init) {
		splitU0NR<LOGN>(U0r, U0i);
    	for (int n = 0; n < LOGN - 1; n++) {
			U0r[n].transpose();
			U0i[n].transpose();
			U0i[n].negate();
			U0r[n] *= (n == 0) ? 0.25 : 0.5;
			U0i[n] *= (n == 0) ? 0.25 : 0.5;
		}
		is_init = true;
	}

	R_Q_square<LOGQ, N> U0ct;
	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA, 3, LOGN-1, G>(U0r, U0i, ct, s, U0ct);

	R_Q_square<LOGQ, N> U0ct_conj;
	int s_conj[N];
	R_Q_square<2*LOGQ, 1 << LOGN> ckey;
	conj<N>(s, s_conj);
	HEAAN<LOGQ,N>::swkgen(s_conj, s, ckey);
	conj(U0ct, ckey, U0ct_conj);

	R_Q<LOGQ,N> pti; pti.setzero();
	pti[N/2][0] = 1; pti[N/2].negate();
	ct_[0] = U0ct; ct_[0] += U0ct_conj;
	ct_[1] = U0ct; ct_[1] -= U0ct_conj; ct_[1]*=pti;
}

//---------------------------------------------
// CoeffToSlot with splitU0 with N=16
//---------------------------------------------
template<int LOGQ>
void CoeffToSlot(const R_Q_square<  LOGQ, 16>& ct,
				 const R_Q_square<2 * LOGQ, 16>& rkey1,
				 const R_Q_square<2 * LOGQ, 16>& rkey2,
				 const R_Q_square<2 * LOGQ, 16>& rkey4,
				 const R_Q_square<2 * LOGQ, 16>& rkey6,
				 const R_Q_square<2 * LOGQ, 16>& rkey7,
				 const R_Q_square<2 * LOGQ, 16>& ckey,
					   R_Q_square<  LOGQ, 16> ct_[2]) {
	SparseDiagonal<8, 3> U0r[3];
	SparseDiagonal<8, 3> U0i[3];
	splitU0NR<4>(U0r, U0i);

	for (int n = 0; n < 3; n++) {
		U0r[n].transpose();
		U0i[n].transpose();
		U0i[n].negate();
	}
	U0r[0] *= (1. / 16);
	U0i[0] *= (1. / 16);

	const R_Q_square<2 * LOGQ, 16>* rkey[3];
	const R_Q_square<    LOGQ, 16> ct1, ct2, temp;
	rkey[0] = &rkey1;
	rkey[1] = &rkey1;
	rkey[2] = &rkey1;
	linear_transform<LOGQ, 16, 50, 3>(U0r[0], U0i[0], ct, rkey[3], ct1);

	rkey[0] = &rkey2;
	rkey[1] = &rkey2;
	rkey[2] = &rkey6;
	linear_transform<LOGQ, 16, 50, 3>(U0r[1], U0i[1], ct1, rkey[3], temp);
		
	rkey[0] = &rkey1;
	rkey[1] = &rkey1;
	rkey[2] = &rkey7;
	linear_transform<LOGQ, 16, 50, 3>(U0r[2], U0i[2], temp, rkey[3], ct1);

	conj<LOGQ, 16>(ct1, ckey, ct2);
	R_Q<LOGQ, 16>pti; pti.setzero();
	pti[8][0] = 1; pti[8].negate();
	ct_[1] = ct1; ct_[1] -= ct2;
	ct_[0] = ct1; ct_[0] -= ct2;

}


template<int LOGQ, int LOGN>
void SlotToCoeff( const R_Q_square<  LOGQ,1<<LOGN>&  ct0,
				  const R_Q_square<  LOGQ,1<<LOGN>&  ct1,
				  const R_Q_square<2*LOGQ,1<<LOGN>& rkey,
					    R_Q_square<  LOGQ,1<<LOGN>& ct_ ){
	const int N = 1 << LOGN;
	double U0r[N/2][N/2];
	double U0i[N/2][N/2]; get_U0<N>(U0r,U0i);
	linear_transform<LOGQ,LOGN,50>(U0r,U0i,ct0,rkey,ct_);

	double iU0r[N/2][N/2];
	double iU0i[N/2][N/2];

	for(int i=0; i<N/2; i++)
	for(int j=0; j<N/2; j++){
		iU0r[i][j] =-U0i[i][j];
		iU0i[i][j] = U0r[i][j];
	}
	R_Q_square<LOGQ,N> temp;
	linear_transform<LOGQ,LOGN,50>(iU0r,iU0i,ct1,rkey,temp);
	ct_+=temp;
}

template<int LOGQ, int LOGN, int LOGDELTA, int G>
void SlotToCoeff( const R_Q_square<  LOGQ,1<<LOGN>&  ct0,
				  const R_Q_square<  LOGQ,1<<LOGN>&  ct1,
				  const int skey[1<<LOGN],
					    R_Q_square<  LOGQ,1<<LOGN>& ct_ ){	
	const int N = 1 << LOGN;
	static bool is_init = false;
	static SparseDiagonal<N/2,3> U0r[LOGN-1];
    static SparseDiagonal<N/2,3> U0i[LOGN-1];
	static SparseDiagonal<N/2,3> iU0r[LOGN-1];
    static SparseDiagonal<N/2,3> iU0i[LOGN-1];

	if(!is_init) {
		SparseDiagonal<N/2,3> U0r_temp[LOGN-1];
    	SparseDiagonal<N/2,3> U0i_temp[LOGN-1];

		splitU0NR<LOGN>(U0r_temp, U0i_temp);
		for(int i = 0; i < LOGN - 1; ++i) {
			U0r[i] = U0r_temp[LOGN - 2 - i];
			U0i[i] = U0i_temp[LOGN - 2 - i];
		}

		for(int i = 0; i < LOGN - 1; ++i) {
			if(i == 0) {
				iU0r[i] = U0i_temp[LOGN - 2 - i];
				iU0r[i].negate();
				iU0i[i] = U0r_temp[LOGN - 2 - i];
			} else {
				iU0r[i] = U0r_temp[LOGN - 2 - i];
				iU0i[i] = U0i_temp[LOGN - 2 - i];
			}
		}		
		is_init = true;
	}

	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA, 3, LOGN-1, G>(U0r, U0i, ct0, skey, ct_);
	R_Q_square<LOGQ, N> ct2;
	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA, 3, LOGN-1, G>(iU0r, iU0i, ct1, skey, ct2);
	ct_ += ct2;
}


template<int LOGN>
void split_U0( double Ar[LOGN][1<<(LOGN-1)][1<<(LOGN-1)],
			   double Ai[LOGN][1<<(LOGN-1)][1<<(LOGN-1)] ){
	int N = 1<<LOGN;
	for(int n=0; n<LOGN; n++)
	for(int i=0; i<N/2;  i++)
	for(int j=0; j<N/2;  j++){ Ar[n][i][j] = 0;
							   Ai[n][i][j] = 0;}
	for(int n=0; n<LOGN-1; n++){
		int M = N>>(n+2);
		for(int i=0, fivei=1<<n; i<M; i++, fivei = (fivei*5)%(2*N)){
			Ar[n][i][i  ] = 1; Ar[n][i+M][i] = 1;
			Ar[n][i][i+M] = cos(PI/N*fivei); Ar[n][i+M][i+M] = -Ar[n][i][i+M];
			Ai[n][i][i+M] = sin(PI/N*fivei); Ai[n][i+M][i+M] = -Ai[n][i][i+M];
		}									   
		for(int k=1; k<(1<<n); k++){ int off = k*2*M;
		for(int i=0; i<M; i++){
			Ar[n][i  +off][i  +off] = Ar[n][i  ][i  ];
			Ar[n][i  +off][i+M+off] = Ar[n][i  ][i+M];
			Ar[n][i+M+off][i  +off] = Ar[n][i+M][i  ];
			Ar[n][i+M+off][i+M+off] = Ar[n][i+M][i+M];
			Ai[n][i  +off][i+M+off] = Ai[n][i  ][i+M];
			Ai[n][i+M+off][i+M+off] = Ai[n][i+M][i+M];
		}}
	}
	for(int i=0; i<N/2; i++){
		int j=0; int i_=i;
		for(int k=0; k<LOGN; k++){
			j+= (i_&1)<<(LOGN-2-k); 
			i_>>=1;
		}

		Ar[LOGN-1][i][j] = 1;
	}
}

template<int LOGQ, int N, int LOGDELTA, int K>
void EvalMod (const R_Q_square<LOGQ, N>& ct, const int s[N], R_Q_square<LOGQ-12*LOGDELTA,N>& p) {
    int s_sq[N];
	conv<N>(s, s, s_sq);

	R_Q_square<2*(LOGQ-   LOGDELTA),N> evk1;
    R_Q_square<2*(LOGQ- 2*LOGDELTA),N> evk2;
    R_Q_square<2*(LOGQ- 3*LOGDELTA),N> evk3;
    R_Q_square<2*(LOGQ- 4*LOGDELTA),N> evk4;
    R_Q_square<2*(LOGQ- 5*LOGDELTA),N> evk5;
    R_Q_square<2*(LOGQ- 6*LOGDELTA),N> evk6;
    R_Q_square<2*(LOGQ- 7*LOGDELTA),N> evk7;
    R_Q_square<2*(LOGQ- 8*LOGDELTA),N> evk8;
    R_Q_square<2*(LOGQ- 9*LOGDELTA),N> evk9;
    R_Q_square<2*(LOGQ-10*LOGDELTA),N> evk10;
    R_Q_square<2*(LOGQ-11*LOGDELTA),N> evk11;

	HEAAN<LOGQ-   LOGDELTA,N>::swkgen(s_sq,s,evk1);
    HEAAN<LOGQ- 2*LOGDELTA,N>::swkgen(s_sq,s,evk2);
    HEAAN<LOGQ- 3*LOGDELTA,N>::swkgen(s_sq,s,evk3);
    HEAAN<LOGQ- 4*LOGDELTA,N>::swkgen(s_sq,s,evk4);
    HEAAN<LOGQ- 5*LOGDELTA,N>::swkgen(s_sq,s,evk5);
    HEAAN<LOGQ- 6*LOGDELTA,N>::swkgen(s_sq,s,evk6);
    HEAAN<LOGQ- 7*LOGDELTA,N>::swkgen(s_sq,s,evk7);
    HEAAN<LOGQ- 8*LOGDELTA,N>::swkgen(s_sq,s,evk8);
    HEAAN<LOGQ- 9*LOGDELTA,N>::swkgen(s_sq,s,evk9);
    HEAAN<LOGQ-10*LOGDELTA,N>::swkgen(s_sq,s,evk10);
    HEAAN<LOGQ-11*LOGDELTA,N>::swkgen(s_sq,s,evk11);

	EvalMod<LOGQ,N,LOGDELTA,K>(ct,evk1,evk2,evk3,evk4,evk5,evk6,evk7,evk8,evk9,evk10,evk11,p);
}

template<int LOGQ, int N, int LOGDELTA, int K>
void EvalMod( const R_Q_square<LOGQ, N>& ct,
		      const R_Q_square<2*(LOGQ-   LOGDELTA), N>& evk1,  // poly
		      const R_Q_square<2*(LOGQ- 2*LOGDELTA), N>& evk2,  // poly
		      const R_Q_square<2*(LOGQ- 3*LOGDELTA), N>& evk3,  // poly
		      const R_Q_square<2*(LOGQ- 4*LOGDELTA), N>& evk4,  // poly
		      const R_Q_square<2*(LOGQ- 5*LOGDELTA), N>& evk5,  // poly
		      const R_Q_square<2*(LOGQ- 6*LOGDELTA), N>& evk6,  // double angle
		      const R_Q_square<2*(LOGQ- 7*LOGDELTA), N>& evk7,  // double angle
		      const R_Q_square<2*(LOGQ- 8*LOGDELTA), N>& evk8,  // double angle
		      const R_Q_square<2*(LOGQ- 9*LOGDELTA), N>& evk9,  // double angle
		      const R_Q_square<2*(LOGQ-10*LOGDELTA), N>& evk10, // arcsine
		      const R_Q_square<2*(LOGQ-11*LOGDELTA), N>& evk11, // arcsine
			        R_Q_square<  (LOGQ-12*LOGDELTA), N>& p    )
{
   	// ct1 = ct * (1/K)
    R_Q_square<LOGQ,N> ct_copy=ct; ct_copy*=uint64_t((1ULL<<LOGDELTA)/static_cast<double>(K));
	R_Q_square<LOGQ-LOGDELTA,N> ct1;
	RS<LOGQ,LOGQ-LOGDELTA,N>(ct_copy, ct1);

    // shift : ct1 -= (1/(4K))	
	Z_Q<LOGQ-LOGDELTA> qtK; qtK.setzero(); qtK[0]=((1ULL<<LOGDELTA)/static_cast<double>(4 * K));
	ct1[0][0] -= qtK;
	
	// Evaluate cos(3 * pi * t) for t in [-1, 1]
    const double u[32] {
-0.18121145350892778,
0.00000000000000176,
-0.51003675488366396,
-0.00000000000000206,
-0.62436800964507644,
-0.00000000000000028,
0.17069002507116071,
0.00000000000000187,
0.64387933774156236,
0.00000000000000155,
-0.60946403046665876,
-0.00000000000000439,
0.16023104261940782,
-0.00000000000000053,
-0.05094267475677791,
0.00000000000000128,
0.00145283705186850,
0.00000000000000268,
-0.00024484541895481,
-0.00000000000000037,
0.00001620811709676,
-0.00000000000000323,
-0.00000171679432516,
0.00000000000000197,
0.00000003753839882,
0.00000000000000263,
-0.00000000274617457,
0.00000000000000504,
0.00000000008565582,
0.00000000000001545,
-0.00000000000460528,
-0.00000000000002388};
	R_Q_square<LOGQ-6*LOGDELTA,N> ct2;
	eval_poly_deg31<LOGQ-LOGDELTA,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,ct2);

	// double angle 1
	R_Q_square<LOGQ-6*LOGDELTA,N> temp1; Mul<LOGQ-6*LOGDELTA,N>(ct2,ct2,evk6,temp1); temp1*=2;
    Z_Q<LOGQ-6*LOGDELTA> one1; one1.setzero(); one1[0]=1ULL<<LOGDELTA; one1*=1ULL<<LOGDELTA;
    temp1[0][0]-=one1;
	R_Q_square<LOGQ-7*LOGDELTA,N> ct3;
	RS<LOGQ-6*LOGDELTA,LOGQ-7*LOGDELTA,N>(temp1, ct3);
	
	// double angle 2
	R_Q_square<LOGQ-7*LOGDELTA,N> temp2; Mul<LOGQ-7*LOGDELTA,N>(ct3,ct3,evk7,temp2); temp2*=2;
    Z_Q<LOGQ-7*LOGDELTA> one2; one2.setzero(); one2[0]=1ULL<<LOGDELTA; one2*=1ULL<<LOGDELTA;
    temp2[0][0]-=one2;
	R_Q_square<LOGQ-8*LOGDELTA,N> ct4;
	RS<LOGQ-7*LOGDELTA,LOGQ-8*LOGDELTA,N>(temp2, ct4);

	// double angle 3
	R_Q_square<LOGQ-8*LOGDELTA,N> temp3; Mul<LOGQ-8*LOGDELTA,N>(ct4,ct4,evk8,temp3); temp3*=2;
    Z_Q<LOGQ-8*LOGDELTA> one3; one3.setzero(); one3[0]=1ULL<<LOGDELTA; one3*=1ULL<<LOGDELTA;
    temp3[0][0]-=one3;
	R_Q_square<LOGQ-9*LOGDELTA,N> ct5;
	RS<LOGQ-8*LOGDELTA,LOGQ-9*LOGDELTA,N>(temp3, ct5);

	// double angle 4
	R_Q_square<LOGQ-9*LOGDELTA,N> temp4; Mul<LOGQ-9*LOGDELTA,N>(ct5,ct5,evk9,temp4); temp4*=2;
    Z_Q<LOGQ-9*LOGDELTA> one4; one4.setzero(); one4[0]=1ULL<<LOGDELTA; one4*=1ULL<<LOGDELTA;
    temp4[0][0]-=one4;
	R_Q_square<LOGQ-10*LOGDELTA,N> ct6;
	RS<LOGQ-9*LOGDELTA,LOGQ-10*LOGDELTA,N>(temp4, ct6);

	// arcsine of degree 3
	// ct7 = x^2 + 6
	R_Q_square<LOGQ-10*LOGDELTA,N> temp5; Mul<LOGQ-10*LOGDELTA,N>(ct6,ct6,evk10,temp5);
	Z_Q<LOGQ-10*LOGDELTA> six; six.setzero(); six[0]=(1ULL<<LOGDELTA)*6; six*=1ULL<<LOGDELTA;
	temp5[0][0]+=six;
	R_Q_square<LOGQ-11*LOGDELTA,N> ct7;
	RS<LOGQ-10*LOGDELTA,LOGQ-11*LOGDELTA,N>(temp5,ct7);
	// ct8 = (1/12pi) * x
	temp5=ct6; temp5 *= static_cast<uint64_t>((1ULL<<LOGDELTA) / (12 * PI));
	R_Q_square<LOGQ-11*LOGDELTA,N> ct8;
	RS<LOGQ-10*LOGDELTA,LOGQ-11*LOGDELTA,N>(temp5,ct8);
	// ct9 = (1/2pi) * (x + (1/6) * x^3)
	R_Q_square<LOGQ-11*LOGDELTA,N> ct9; 
	Mul<LOGQ-11*LOGDELTA,N>(ct7,ct8,evk11,ct9);
	RS<LOGQ-11*LOGDELTA,LOGQ-12*LOGDELTA,N>(ct9,p);	    
}

// Rescale by the size of q should not be performed at the beginning of EvalMod_Kx.
/*
template<int LOGQ, int N, int LOGq, int LOGDELTA>
void EvalMod_K4( const R_Q_square<  LOGQ,N>& ct,
				  const R_Q_square<2*(LOGQ-LOGq-2         ),N>& evk1,
				  const R_Q_square<2*(LOGQ-LOGq-2-  LOGDELTA),N>& evk2,
				  const R_Q_square<2*(LOGQ-LOGq-2-2*LOGDELTA),N>& evk3,
				  const R_Q_square<2*(LOGQ-LOGq-2-3*LOGDELTA),N>& evk4,
			      const R_Q_square<2*(LOGQ-LOGq-2-4*LOGDELTA),N>& evk5,
						R_Q_square<  (LOGQ-LOGq-2-5*LOGDELTA),N>& p    )
{	R_Q_square<LOGQ-LOGq-2,N> ct1;
	RS<LOGQ,LOGQ-LOGq-2,N>(ct,ct1); // Why rescale?
	ct1.print();
	const double u[32] = { 0.000000000000,0.140829338704,0.000000000000,-0.071956680803,
						   0.000000000000,0.314734688318,0.000000000000,-0.614678201819,
					      -0.000000000000,0.174769865184,-0.000000000000,0.037343013007,
					      -0.000000000000,-0.083190952031,-0.000000000000,0.128405049045,
						  -0.000000000000,0.247970475294,0.000000000000,-0.195489024287,
						   0.000000000000,0.784513028614,0.000000000000,-0.995925646437,
						   0.000000000000,0.336175543041,0.000000000000,-0.246866274711,
						 -0.000000000000,0.083032826731,-0.000000000000,-0.040430268293};
	eval_poly_deg31<LOGQ-LOGq-2,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,p);
	p.print();
	p *= 1ULL<<(LOGq+2);
}


template<int LOGQ, int N, int LOGq, int LOGDELTA>
void EvalMod_K2( const R_Q_square<  LOGQ,N>& ct,
				  const R_Q_square<2*(LOGQ-LOGq-1         ),N>& evk1,
				  const R_Q_square<2*(LOGQ-LOGq-1-  LOGDELTA),N>& evk2,
				  const R_Q_square<2*(LOGQ-LOGq-1-2*LOGDELTA),N>& evk3,
				  const R_Q_square<2*(LOGQ-LOGq-1-3*LOGDELTA),N>& evk4,
			      const R_Q_square<2*(LOGQ-LOGq-1-4*LOGDELTA),N>& evk5,
						R_Q_square<  (LOGQ-LOGq-1-5*LOGDELTA),N>& p    )
{	R_Q_square<LOGQ-LOGq-1,N> ct1;
	RS<LOGQ,LOGQ-LOGq-1,N>(ct,ct1);
	ct1.print();
	const double u[32] = {0.000000000000,0.132748607764,0.000000000000,-0.380278105603,0.000000000000,0.162806166803,-0.000000000000,0.107579389804,-0.000000000000,0.417564594173,-0.000000000000,-0.575577168110,-0.000000000000,0.279185145530,-0.000000000000,-0.149544709365,-0.000000000000,0.007450388247,0.000000000000,-0.002135232830,0.000000000000,0.000243436777,-0.000000000000,-0.000043966660,0.000000000000,0.000001649668,0.000000000000,-0.000000206220,-0.000000000000,0.000000011034,-0.000000000000,-0.000000001015 };
	eval_poly_deg31<LOGQ-LOGq-1,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,p);
	p.print();
	p *= 1ULL<<(LOGq);
}

template<int LOGQ, int N, int LOGq, int LOGDELTA>
void EvalMod_K3( const R_Q_square<  LOGQ,N>& ct,
				  const R_Q_square<2*(LOGQ-LOGq-  LOGDELTA),N>& evk1,
				  const R_Q_square<2*(LOGQ-LOGq-2*LOGDELTA),N>& evk2,
				  const R_Q_square<2*(LOGQ-LOGq-3*LOGDELTA),N>& evk3,
				  const R_Q_square<2*(LOGQ-LOGq-4*LOGDELTA),N>& evk4,
			      const R_Q_square<2*(LOGQ-LOGq-5*LOGDELTA),N>& evk5,
						R_Q_square<  (LOGQ-LOGq-6*LOGDELTA),N>& p    )
{	
	R_Q_square<  LOGQ,N> ct_copy=ct; ct_copy*=uint64_t((1ULL<<LOGDELTA)/3.);
	R_Q_square<LOGQ-LOGq-LOGDELTA,N> ct1;
	RS<LOGQ,LOGQ-LOGq-LOGDELTA,N>(ct_copy,ct1);
	ct1.print();
	const double u[32] = {-0.000000000000,0.137852010310,0.000000000000,-0.018904701355,0.000000000000,0.024320661898,-0.000000000000,-0.151590145511,-0.000000000000,0.367515735123,-0.000000000000,0.056870745879,-0.000000000000,0.541162271614,-0.000000000000,-1.048431384285,-0.000000000000,0.316696874787,0.000000000000,-0.282534155918,0.000000000000,0.106543202085,-0.000000000000,-0.053509305951,0.000000000000,0.005651171532,0.000000000000,-0.001841556587,-0.000000000000,0.000257314691,-0.000000000000,-0.000059365641};
	eval_poly_deg31<LOGQ-LOGq-LOGDELTA,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,p);
	p.print();
	p *= 1ULL<<(LOGq);
}
*/
