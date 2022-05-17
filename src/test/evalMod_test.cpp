#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"
#include "util/util.h"

#include "experiment/rns_debug.h"

#include <iostream>

void EvalMod_test()
{
    int s[N], s_sq[N];
    HEAAN<LOGQ,N>::keygen(H,s);
	conv<N>(s, s, s_sq);

	R_Q_square<2*(LOGQ-   LOGDELTA_boot),N> evk1;
	R_Q_square<2*(LOGQ- 2*LOGDELTA_boot),N> evk2;
	R_Q_square<2*(LOGQ- 3*LOGDELTA_boot),N> evk3;
	R_Q_square<2*(LOGQ- 4*LOGDELTA_boot),N> evk4;
	R_Q_square<2*(LOGQ- 5*LOGDELTA_boot),N> evk5;
	R_Q_square<2*(LOGQ- 6*LOGDELTA_boot),N> evk6;
	R_Q_square<2*(LOGQ- 7*LOGDELTA_boot),N> evk7;
	R_Q_square<2*(LOGQ- 8*LOGDELTA_boot),N> evk8;
	R_Q_square<2*(LOGQ- 9*LOGDELTA_boot),N> evk9;
	R_Q_square<2*(LOGQ-10*LOGDELTA_boot),N> evk10;
	R_Q_square<2*(LOGQ-11*LOGDELTA_boot),N> evk11;

	HEAAN<LOGQ-   LOGDELTA_boot,N>::swkgen(s_sq,s,evk1);
	HEAAN<LOGQ- 2*LOGDELTA_boot,N>::swkgen(s_sq,s,evk2);
	HEAAN<LOGQ- 3*LOGDELTA_boot,N>::swkgen(s_sq,s,evk3);
	HEAAN<LOGQ- 4*LOGDELTA_boot,N>::swkgen(s_sq,s,evk4);
	HEAAN<LOGQ- 5*LOGDELTA_boot,N>::swkgen(s_sq,s,evk5);
	HEAAN<LOGQ- 6*LOGDELTA_boot,N>::swkgen(s_sq,s,evk6);
	HEAAN<LOGQ- 7*LOGDELTA_boot,N>::swkgen(s_sq,s,evk7);
	HEAAN<LOGQ- 8*LOGDELTA_boot,N>::swkgen(s_sq,s,evk8);
	HEAAN<LOGQ- 9*LOGDELTA_boot,N>::swkgen(s_sq,s,evk9);
	HEAAN<LOGQ-10*LOGDELTA_boot,N>::swkgen(s_sq,s,evk10);
	HEAAN<LOGQ-11*LOGDELTA_boot,N>::swkgen(s_sq,s,evk11);

    Message<LOGN> z, z_out;
    set_evalmod_message(z);

    R_Q<LOGQ, N> pt;
    R_Q_square<LOGQ,N> ct;
    encode(z,Delta_boot,pt);
    HEAAN<LOGQ,N>::enc(pt,s,ct);

	R_Q_square<LOGQ-12*LOGDELTA_boot,N> ct_out;

	EvalMod<LOGQ,N,LOGDELTA_boot,K>(ct,evk1,evk2,evk3,evk4,evk5,evk6,evk7,evk8,evk9,evk10,evk11,ct_out);

	R_Q<LOGQ-12*LOGDELTA_boot,N> pt_out;
    HEAAN<LOGQ-12*LOGDELTA_boot,N>::dec(ct_out,s,pt_out);
	decode_log(pt_out, LOGDELTA, z_out);
	print("z", z);
	print("z_out", z_out);
}

/*
void CoeffToSlot_SlotToCoeff_test()
{
    int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);
    
    Message<LOGN> z, z_out;
	set_random_message(z);

    R_Q<LOGQ, N> pt, pt_out;
    R_Q_square<LOGQ,N> ct, ct_cts[2], ct_out;
    encode(z,Delta,pt);
	HEAAN<LOGQ,N>::enc(pt,s,ct);
	
	CoeffToSlot<LOGQ,LOGN, LOGDELTA_boot_tilde, G>(ct,s,ct_cts);
    SlotToCoeff<LOGQ,LOGN, LOGDELTA_boot, G>(ct_cts[0], ct_cts[1],s,ct_out);
    HEAAN<LOGQ,N>::dec(ct_out,s,pt_out);
    decode_log(pt_out,LOGDELTA +(LOGN-1)/G*LOGDELTA_boot_tilde + (LOGN-1)/G*LOGDELTA_boot,z_out);
    print("z", z);
    print("z_out", z_out);
}
*/

int main()
{
	EvalMod_test();
}
