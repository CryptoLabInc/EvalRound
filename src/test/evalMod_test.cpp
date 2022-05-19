#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"
#include "util/util.h"

#include "experiment/rns_debug.h"

#include <iostream>

void EvalMod_test()
{
    int s[N], s_sq[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z, z_out;
    set_evalmod_message(z);

    R_Q<LOGQ, N> pt;
    R_Q_square<LOGQ,N> ct;
    encode(z,Delta_boot,pt);
    HEAAN<LOGQ,N>::enc(pt,s,ct);

	R_Q_square<LOGQ-12*LOGDELTA_boot,N> ct_out;

	EvalMod<LOGQ,N,LOGDELTA_boot,K>(ct,s,ct_out);

	R_Q<LOGQ-12*LOGDELTA_boot,N> pt_out;
    HEAAN<LOGQ-12*LOGDELTA_boot,N>::dec(ct_out,s,pt_out);
	decode_log(pt_out, LOGDELTA_boot, z_out);
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
