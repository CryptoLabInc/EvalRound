#include "HEAAN/endecode.h"
#include "HEAAN/arith/conv.h"
#include "experiment/test.h"

template<int LOGN>
void ctxt_conv_test(){	
    const int N = 1 << LOGN;
	Message<LOGN> z, z_sq_orig, z_sq;
	set_random_message(z);
    mul(z, z, z_sq_orig);

	int s[N], s_sq[N];
    R_Q_square<2*LOGQ,N>evk;
    HEAAN<LOGQ,N>::keygen(H,s);
    conv<N>(s, s, s_sq);
	HEAAN<LOGQ,N>::swkgen(s_sq,s,evk);

    R_Q<LOGQ, N> pt, pt_sq;
    R_Q_square<LOGQ,N> ct, ct_sq;
    encode(z,Delta,pt);
	HEAAN<LOGQ,N>::enc(pt,s,ct);
	Mul(ct,ct,evk,ct_sq);
	HEAAN<LOGQ,N>::dec(ct_sq,s,pt_sq);
	decode_log(pt_sq,2*LOGDELTA,z_sq);

    print("z_sq_orig", z_sq_orig);
    print("z_sq", z_sq);
}

int main() {
    ctxt_conv_test<9>();
}