#pragma once

#include <cmath>

#include "linear_transform.h"
#include "poly.h"

template< int LOGQfr, int LOGQto, int N >
void mod_raise(const R_Q_square<LOGQfr,N>& ctfr,
					 R_Q_square<LOGQto,N>& ctto){
	resize(ctfr, ctto);
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
			U0r[n] *= 1 / pow(N, 1.0 / (LOGN - 1));
			U0i[n] *= 1 / pow(N, 1.0 / (LOGN - 1));
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