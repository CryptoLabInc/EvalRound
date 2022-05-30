#include "HEAAN/bootstrap.h"
#include "experiment/rns_debug.h"

#include <iostream>

void bootstrap_test()
{
    int s[N], s_sq[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z, z_out;
    set_test_message(z);

    R_Q<LOGq, N> pt;
    R_Q_square<LOGq,N> ct;
    encode(z,Delta,pt);
    HEAAN<LOGq,N>::enc(pt,s,ct);

	// ModRaise
	R_Q_square<LOGQ,N> ct_modraise;
	mod_raise<LOGq,LOGQ,N>(ct,ct_modraise);

	// CoeffToSlot
	R_Q_square<LOGQ,N> ct_cts[2];
	CoeffToSlot<LOGQ,LOGN,LOGDELTA_boot_tilde,G>(ct_modraise,s,ct_cts);
	const int LOGQ_after_cts = LOGQ-(LOGN-1)/G*LOGDELTA_boot_tilde;
	R_Q_square<LOGQ_after_cts,N> ct_ctsrs[2];
	for(int i=0;i<2;i++)
	    RS<LOGQ,LOGQ_after_cts,N>(ct_cts[i],ct_ctsrs[i]);

	// EvalMod
	const int LOGQ_after_evalmod = LOGQ_after_cts - 12*LOGDELTA_boot;
	R_Q_square<LOGQ_after_evalmod,N> ct_evalmod[2];
	for(int i=0;i<2;i++)
	    EvalMod<LOGQ_after_cts,N,LOGDELTA_boot,K>(ct_ctsrs[i],s,ct_evalmod[i]);

	// SlotToCoeff
	// STC is done on LOGDELTA since it doesn't have to handle (pt + qI)
	const int LOGQ_after_stc = LOGQ_after_evalmod - (LOGN-1)/G*LOGDELTA;
    R_Q_square<LOGQ_after_evalmod,N> ct_stc;
	R_Q_square<LOGQ_after_stc,N> ct_boot;
	SlotToCoeff<LOGQ_after_evalmod,LOGN,LOGDELTA,G>(ct_evalmod[0],ct_evalmod[1],s,ct_stc);
	RS<LOGQ_after_evalmod,LOGQ_after_stc,N>(ct_stc,ct_boot);

	R_Q<LOGQ_after_stc,N> pt_out;
	HEAAN<LOGQ_after_stc,N>::dec(ct_boot,s,pt_out);
	decode_log(pt_out,LOGDELTA,z_out);
	print("z", z);
	print("z_out", z_out);
	print_max_error<LOGN>(z,z_out);

	R_Q<LOGQ_after_stc, N> e;
	resize(pt, e);
	e -= pt_out;

	double e_per_Delta_sup_norm = 0;
    for(int i = 0; i < N; ++i) {
        Z_Q<LOGQ_after_stc> val = e[i];
		if(val.is_bigger_than_halfQ())
            val.negate();
        double val_abs_double = (double) (val[0]) / Delta;

        e_per_Delta_sup_norm = val_abs_double > e_per_Delta_sup_norm ? val_abs_double : e_per_Delta_sup_norm;
    }
	std::cout << "LOG2 sup_norm(pt/Delta - pt_tilde/Delta) : " << std::log2(e_per_Delta_sup_norm) << std::endl;
}

int main()
{
	bootstrap_test();
}
