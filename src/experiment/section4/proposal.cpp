#include "HEAAN/bootstrap.h"
#include "experiment/rns.h"

#include <iostream>

template<int LOGN, int LOGDELTA_cts, int LOGDELTA_stc, int G>
void conventional_bootstrap_test()
{
    std::cout << "Measuring error on conventional bootstrap" << std::endl;
	std::cout << "LOGN : " << LOGN << std::endl;
	std::cout << "LOGDELTA_cts : " << LOGDELTA_cts << std::endl;
	std::cout << "LOGDELTA_stc : " << LOGDELTA_stc << std::endl;

	const int N = 1 << LOGN;

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
	CoeffToSlot<LOGQ,LOGN,LOGDELTA_cts,G>(ct_modraise,s,ct_cts);
	const int LOGQ_after_cts = LOGQ-(LOGN-1)/G*LOGDELTA_cts;
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
	const int LOGQ_after_stc = LOGQ_after_evalmod - (LOGN-1)/G*LOGDELTA_stc;
    R_Q_square<LOGQ_after_evalmod,N> ct_stc;
	R_Q_square<LOGQ_after_stc,N> ct_boot;
	SlotToCoeff<LOGQ_after_evalmod,LOGN,LOGDELTA_stc,G>(ct_evalmod[0],ct_evalmod[1],s,ct_stc);
	RS<LOGQ_after_evalmod,LOGQ_after_stc,N>(ct_stc,ct_boot);

	R_Q<LOGQ_after_stc,N> pt_out;
	HEAAN<LOGQ_after_stc,N>::dec(ct_boot,s,pt_out);
	decode_log(pt_out,LOGDELTA,z_out);
	print_max_error<LOGN>(z,z_out);

	R_Q<LOGQ_after_stc, N> e;
	resize(pt, e);
	e -= pt_out;

	double e_per_Delta_sup_norm = 0;
    for(int i = 0; i < N; ++i) {
		double val = (double) e[i];
		double val_abs_double = (double) std::abs(val) / Delta;

        e_per_Delta_sup_norm = val_abs_double > e_per_Delta_sup_norm ? val_abs_double : e_per_Delta_sup_norm;
    }
	std::cout << "sup_norm(pt/Delta - pt_tilde/Delta) : " << e_per_Delta_sup_norm << std::endl;
	std::cout << "Modulus consumed : " << (LOGQ - LOGQ_after_stc) << std::endl;
}

template<int LOGN, int LOGDELTA_cts, int LOGDELTA_stc, int G>
void proposed_bootstrap_test()
{
    std::cout << "Measuring error on proposed bootstrap" << std::endl;
	std::cout << "LOGN : " << LOGN << std::endl;
	std::cout << "LOGDELTA_cts : " << LOGDELTA_cts << std::endl;
	std::cout << "LOGDELTA_stc : " << LOGDELTA_stc << std::endl;

    const int N = 1 << LOGN;

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
	CoeffToSlot<LOGQ,LOGN,LOGDELTA_cts,G>(ct_modraise,s,ct_cts);
	const int LOGQ_after_cts = LOGQ-(LOGN-1)/G*LOGDELTA_cts;
	R_Q_square<LOGQ_after_cts,N> ct_ctsrs[2];
	for(int i=0;i<2;i++)
	    RS<LOGQ,LOGQ_after_cts,N>(ct_cts[i],ct_ctsrs[i]);

	// EvalMod
	const int LOGQ_after_evalmod = LOGQ_after_cts - 12*LOGDELTA_boot;
	R_Q_square<LOGQ_after_evalmod,N> ct_evalmod[2];
	for(int i=0;i<2;i++)
	    EvalMod<LOGQ_after_cts,N,LOGDELTA_boot,K>(ct_ctsrs[i],s,ct_evalmod[i]);
	
	// EvalqI
	R_Q_square<LOGQ_after_evalmod,N> ct_evalqI[2];
	for(int i = 0; i < 2; i++) {
		resize(ct_ctsrs[i], ct_evalqI[i]);
		ct_evalqI[i] -= ct_evalmod[i];
	}

	// SlotToCoeff
	const int LOGQ_after_stc = LOGQ_after_evalmod - (LOGN-1)/G*LOGDELTA_stc;
    R_Q_square<LOGQ_after_evalmod,N> ct_stc;
	R_Q_square<LOGQ_after_stc,N> ct_qI;
	SlotToCoeff<LOGQ_after_evalmod,LOGN,LOGDELTA_stc,G>(ct_evalqI[0],ct_evalqI[1],s,ct_stc);
	RS<LOGQ_after_evalmod,LOGQ_after_stc,N>(ct_stc, ct_qI);
	
	// Sub
	R_Q_square<LOGQ_after_stc,N> ct_boot;
	resize(ct_modraise, ct_boot);
	ct_boot -= ct_qI;

	R_Q<LOGQ_after_stc,N> pt_out;
	HEAAN<LOGQ_after_stc,N>::dec(ct_boot,s,pt_out);
	decode_log(pt_out,LOGDELTA,z_out);
	print_max_error<LOGN>(z,z_out);

	R_Q<LOGQ_after_stc, N> e;
	resize(pt, e);
	e -= pt_out;

	double e_per_Delta_sup_norm = 0;
    for(int i = 0; i < N; ++i) {
    	double val = (double) e[i];
		double val_abs_double = (double) std::abs(val) / Delta;

        e_per_Delta_sup_norm = val_abs_double > e_per_Delta_sup_norm ? val_abs_double : e_per_Delta_sup_norm;
    }
	std::cout << "sup_norm(pt/Delta - pt_tilde/Delta) : " << e_per_Delta_sup_norm << std::endl;
	std::cout << "Modulus consumed : " << (LOGQ - LOGQ_after_stc) << std::endl;
}

int main()
{
    conventional_bootstrap_test<9, 22, 50, 2>(); proposed_bootstrap_test<9, 22, 60, 2>();
	conventional_bootstrap_test<9, 60, 50, 2>(); proposed_bootstrap_test<9, 60, 60, 2>();
	conventional_bootstrap_test<13, 25, 50, 3>(); proposed_bootstrap_test<13, 25, 60, 3>();
	conventional_bootstrap_test<13, 60, 50, 3>(); proposed_bootstrap_test<13, 60, 60, 3>();
    conventional_bootstrap_test<17, 29, 50, 4>(); proposed_bootstrap_test<17, 29, 60, 4>();
	conventional_bootstrap_test<17, 60, 50, 4>(); proposed_bootstrap_test<17, 60, 60, 4>();
}