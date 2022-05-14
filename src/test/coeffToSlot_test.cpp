#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"
#include "util/util.h"

#include "experiment/big.h"

#include <iostream>

int main()
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
    //SlotToCoeff<LOGQ,LOGN, LOGDELTA>(ct_cts[0], ct_[1],s,ct);
    decode_log(pt_cts,LOGDELTA +9*LOGDELTA_TILDE,z_cts);

    for(int i = 0; i < 10; ++i) {
        Z_Q<LOGQ> val(pt[bitReverse(i, LOGN - 1)]);
        bool is_negative =  val.is_bigger_than_halfQ();
        if(is_negative)
            val.negate();
        double val_double = (double) val[0];
        if(is_negative)
            val_double *= -1;   
        val_double /= 2 * Delta;
        //for(int j = 0; j < 3; ++j)
        //std::cout << pt[bitReverse(i, LOGN - 1)][j] << std::endl;
        std::cout << z_cts.r[i] << " " << val_double << std::endl;
    }

}