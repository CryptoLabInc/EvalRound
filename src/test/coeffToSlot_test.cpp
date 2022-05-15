#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"
#include "util/util.h"

#include "experiment/big.h"

#include <iostream>

void CoeffToSlot_test()
{
    int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);
    
    Message<LOGN> z, z_cts;
	set_random_message(z);

    R_Q<LOGQ, N> pt, pt_cts;
    R_Q_square<LOGQ,N> ct, ct_cts[2];
    encode(z,Delta,pt);
	HEAAN<LOGQ,N>::enc(pt,s,ct);
	
	CoeffToSlot<LOGQ,LOGN, LOGDELTA_TILDE>(ct,s,ct_cts);
    HEAAN<LOGQ,N>::dec(ct_cts[0],s,pt_cts);
    const int D = LOGN == 10? 4 : LOGN-1;
    decode_log(pt_cts,LOGDELTA +D*LOGDELTA_TILDE,z_cts);

    for(int i = 0; i < 10; ++i) {
        Z_Q<LOGQ> val(pt[bitReverse(i, LOGN - 1)]);
        bool is_negative =  val.is_bigger_than_halfQ();
        if(is_negative)
            val.negate();
        double val_double = (double) val[0];
        if(is_negative)
            val_double *= -1;   
        val_double /= 2 * Delta;
        std::cout << z_cts.r[i] << " " << val_double << std::endl;
    }
}

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
	
	CoeffToSlot<LOGQ,LOGN, LOGDELTA_TILDE>(ct,s,ct_cts);
    SlotToCoeff<LOGQ,LOGN, LOGDELTA_TILDE>(ct_cts[0], ct_cts[1],s,ct_out);
    HEAAN<LOGQ,N>::dec(ct_out,s,pt_out);
    const int D = LOGN == 10? 4 : LOGN-1;
    decode_log(pt_out,LOGDELTA +2*D*LOGDELTA_TILDE,z_out);
    print("z", z);
    print("z_out", z_out);
}

int main()
{
    CoeffToSlot_SlotToCoeff_test();
}