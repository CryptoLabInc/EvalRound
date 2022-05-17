#include "experiment/rns_debug.h"

#include <iostream>

int main()
{
    int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);
    
    Message<LOGN> z;
	set_random_message(z);

    R_Q<LOGQ, N> pt;
    R_Q_square<LOGQ,N> ct;
    R_Q_square<LOGq,N> ct_q;
    R_Q<LOGq, N> pt_q;
    Message<LOGN> z_q;
    encode(z,Delta,pt);
	HEAAN<LOGQ,N>::enc(pt,s,ct);
    resize(ct, ct_q); // size down to RQ -> Rq
    HEAAN<LOGq,N>::dec(ct_q,s, pt_q);
    decode(pt_q, Delta, z_q);
    print("z", z);
    print("z_q", z_q); // check z_q is not damaged
}