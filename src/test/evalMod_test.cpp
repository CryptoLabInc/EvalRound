#include "HEAAN/bootstrap.h"
#include "experiment/test.h"

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

int main()
{
	EvalMod_test();
}
