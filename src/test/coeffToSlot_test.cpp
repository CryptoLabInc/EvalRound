#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"
#include "util/util.h"

#include "experiment/big.h"

#include <iostream>

void CoeffToSlot_test()
{
    int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);
   
    Message<LOGN> z, z_cts[2], z_cts_res, z_cts_exact, e;
    set_random_message(z);

    R_Q<LOGQ, N> pt, pt_cts[2];
    R_Q_square<LOGQ,N> ct, ct_cts[2];
    encode(z,Delta,pt);
    HEAAN<LOGQ,N>::enc(pt,s,ct);

    CoeffToSlot<LOGQ,LOGN, LOGDELTA_TILDE>(ct,s,ct_cts);
    for(int i = 0; i < 2; ++i) {
        HEAAN<LOGQ,N>::dec(ct_cts[i],s,pt_cts[i]);
        decode_log(pt_cts[i],LOGDELTA +(LOGN-1)/2*LOGDELTA_TILDE,z_cts[i]);
    }
    for(int i = 0; i < N/2; ++i) {
        z_cts_res.r[i] = z_cts[0].r[i];
        z_cts_res.i[i] = z_cts[1].r[i];
    }

    for(int i = 0; i < N/2; ++i) {
        Z_Q<LOGQ> valr(pt[bitReverse(i, LOGN - 1)]);
        Z_Q<LOGQ> vali(pt[bitReverse(i, LOGN - 1) + N/2]);
        bool is_negative_r =  valr.is_bigger_than_halfQ();
        if(is_negative_r)
            valr.negate();
        double valr_double = (double) valr[0];
        if(is_negative_r)
            valr_double *= -1;  
        valr_double /= Delta;

        bool is_negative_i =  vali.is_bigger_than_halfQ();
        if(is_negative_i)
            vali.negate();
        double vali_double = (double) vali[0];
        if(is_negative_i)
            vali_double *= -1;  
        vali_double /= Delta;

        z_cts_exact.r[i] = valr_double;
        z_cts_exact.i[i] = vali_double;
    }


    sub(z_cts_res, z_cts_exact, e);
    print("z_cts_res", z_cts_res);
    print("z_cts_exact", z_cts_exact);
    print("e", e);
    std::cout << norm(e) << " " << norm(z) << " " << (norm(e) / norm(z)) << std::endl;
    std::cout << sup_norm(z_cts_res) << " " << sup_norm(e) << std::endl;
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
    decode_log(pt_out,LOGDELTA +(LOGN-1)*LOGDELTA_TILDE,z_out);
    print("z", z);
    print("z_out", z_out);
}

int main()
{
    CoeffToSlot_test();
    CoeffToSlot_SlotToCoeff_test();
}