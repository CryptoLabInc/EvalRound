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
    encode(z,Delta,pt);
	HEAAN<LOGQ,N>::enc(pt,s,ct);
}