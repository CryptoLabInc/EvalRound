#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"

#include "experiment/big.h"

#include <iostream>

int main()
{
    int s[N], s_rot[N], s_conj[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z, z_measured;
	set_random_message(z);

    R_Q<LOGQ, N> pt;
    R_Q_square<LOGQ,N> ct, ct_[2];
    encode(z,Delta,pt);
	HEAAN<LOGQ,N>::enc(pt,s,ct);
	
	CoeffToSlot<LOGQ,LOGN, LOGDELTA>(ct,s,ct_);
	//SlotToCoeff<LOGQ,LOGN>(ct_[0], ct_[1],rkey,ct);

	HEAAN<LOGQ,N>::dec(ct,s,pt);
	decode_log(pt,LOGDELTA,z_measured);
    print("z", z);
    print("z_measured", z_measured);
}