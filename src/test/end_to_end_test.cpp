#include "HEAAN/bootstrap.h"
#include "experiment/rns.h"

#include <iostream>

template<int LOGN, int LOGDELTA_cts, int LOGDELTA_stc, int G>
void end_to_end_test()
{
	std::cout << LOGN << " " << LOGDELTA_cts << " " << LOGDELTA_stc << " " << G << std::endl;
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

	// CoeffToSlot Test
	{
		Message<LOGN> z_cts[2], z_cts_exact[2], e[2];
		R_Q<LOGQ_after_cts, N> pt_cts[2];
		for(int i = 0; i < 2; ++i) {
        	HEAAN<LOGQ_after_cts,N>::dec(ct_ctsrs[i],s,pt_cts[i]);
        	decode_log(pt_cts[i],LOGDELTA,z_cts[i]);
    	}
	
		// compute desired z_cts
		R_Q<LOGQ, N> pt_modraise;
		HEAAN<LOGQ,N>::dec(ct_modraise,s, pt_modraise);
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt_modraise[bitReverse(i, LOGN - 1)];
			z_cts_exact[0].r[i] = val_double / Delta;
		}
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt_modraise[bitReverse(i, LOGN - 1) + N/2];
			z_cts_exact[1].r[i] = val_double / Delta;
		}

		// get error
	    for(int i = 0; i < 2; ++i)
        	sub(z_cts[i], z_cts_exact[i], e[i]);

		std::cout << "cts error : " << std::max(sup_norm(e[0]), sup_norm(e[1])) << std::endl;
	}

	// EvalMod
	const int LOGQ_after_evalmod = LOGQ_after_cts - 12*LOGDELTA_boot;
	R_Q_square<LOGQ_after_evalmod,N> ct_evalmod[2];
	for(int i=0;i<2;i++)
	    EvalMod<LOGQ_after_cts,N,LOGDELTA_boot,K>(ct_ctsrs[i],s,ct_evalmod[i]);
	
	// EvalMod Test
	{
		Message<LOGN> z_evalmod[2], z_evalmod_exact[2], e[2];
		R_Q<LOGQ_after_evalmod, N> pt_evalmod[2];
		for(int i = 0; i < 2; ++i) {
        	HEAAN<LOGQ_after_evalmod,N>::dec(ct_evalmod[i],s,pt_evalmod[i]);
        	decode_log(pt_evalmod[i],LOGDELTA,z_evalmod[i]);
    	}
	
		// compute desired z_evalmod
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt[bitReverse(i, LOGN - 1)];
			z_evalmod_exact[0].r[i] = val_double / Delta;
		}
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt[bitReverse(i, LOGN - 1) + N/2];
			z_evalmod_exact[1].r[i] = val_double / Delta;
		}

		// get error
	    for(int i = 0; i < 2; ++i)
        	sub(z_evalmod[i], z_evalmod_exact[i], e[i]);

		std::cout << "evalmod error : " << std::max(sup_norm(e[0]), sup_norm(e[1])) << std::endl;
	}

	// EvalqI
	R_Q_square<LOGQ_after_evalmod,N> ct_evalqI[2];
	for(int i = 0; i < 2; i++) {
		resize(ct_ctsrs[i], ct_evalqI[i]);
		ct_evalqI[i] -= ct_evalmod[i];
	}

	// EvalqI Test
	{
		Message<LOGN> z_evalqI[2], z_evalqI_exact[2], e[2];
		R_Q<LOGQ_after_evalmod, N> pt_evalqI[2];
		for(int i = 0; i < 2; ++i) {
        	HEAAN<LOGQ_after_evalmod,N>::dec(ct_evalqI[i],s,pt_evalqI[i]);
        	decode_log(pt_evalqI[i],LOGDELTA,z_evalqI[i]);
    	}
	
		// compute desired z_evalqI
		R_Q<LOGQ, N> pt_modraise;
		HEAAN<LOGQ,N>::dec(ct_modraise,s, pt_modraise);
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt_modraise[bitReverse(i, LOGN - 1)] - (double) pt[bitReverse(i, LOGN - 1)];
			z_evalqI_exact[0].r[i] = val_double / Delta;
		}
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt_modraise[bitReverse(i, LOGN - 1) + N/2] - (double) pt[bitReverse(i, LOGN - 1) + N/2];
			z_evalqI_exact[1].r[i] = val_double / Delta;
		}

		// get error
	    for(int i = 0; i < 2; ++i)
        	sub(z_evalqI[i], z_evalqI_exact[i], e[i]);

		std::cout << "evalqI error : " << std::max(sup_norm(e[0]), sup_norm(e[1])) << std::endl;
	}
	
	// SlotToCoeff
	const int LOGQ_after_stc = LOGQ_after_evalmod - (LOGN-1)/G*LOGDELTA_stc;
    R_Q_square<LOGQ_after_evalmod,N> ct_stc;
	R_Q_square<LOGQ_after_stc,N> ct_qI;
	SlotToCoeff<LOGQ_after_evalmod,LOGN,LOGDELTA_stc,G>(ct_evalqI[0],ct_evalqI[1],s,ct_stc);
	RS<LOGQ_after_evalmod,LOGQ_after_stc,N>(ct_stc, ct_qI);

	// SlotToCoeff Test
	{
		// compute qI by evalqI
		Message<LOGN> z_evalqI[2], z_evalqI_stc[2], e[2];
		R_Q<LOGQ_after_evalmod, N> pt_evalqI[2];
		for(int i = 0; i < 2; ++i) {
        	HEAAN<LOGQ_after_evalmod,N>::dec(ct_evalqI[i],s,pt_evalqI[i]);
        	decode_log(pt_evalqI[i],LOGDELTA,z_evalqI[i]);
    	}

		// after stc
		R_Q<LOGQ_after_stc, N> pt_qI;
		HEAAN<LOGQ_after_stc, N>::dec(ct_qI, s, pt_qI);
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt_qI[bitReverse(i, LOGN - 1)];
			z_evalqI_stc[0].r[i] = val_double / Delta;
		}
		for(int i = 0; i < N/2; ++i) {
			double val_double = (double) pt_qI[bitReverse(i, LOGN - 1) + N/2];
			z_evalqI_stc[1].r[i] = val_double / Delta;
		}

		// get error
	    for(int i = 0; i < 2; ++i)
        	sub(z_evalqI[i], z_evalqI_stc[i], e[i]);

		std::cout << "stc error : " << std::max(sup_norm(e[0]), sup_norm(e[1])) << std::endl;
	}
	
	// Sub
	R_Q_square<LOGQ_after_stc,N> ct_boot;
	resize(ct_modraise, ct_boot);
	ct_boot -= ct_qI;

	R_Q<LOGQ_after_stc,N> pt_out;
	HEAAN<LOGQ_after_stc,N>::dec(ct_boot,s,pt_out);
	decode_log(pt_out,LOGDELTA,z_out);

	{
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
	}

	print_max_error<LOGN>(z,z_out);
	std::cout << std::endl;
}

int main()
{
    end_to_end_test<9, 22, 50, 2>(); end_to_end_test<9, 22, 60, 2>();
	end_to_end_test<9, 60, 50, 2>(); end_to_end_test<9, 60, 60, 2>();
	end_to_end_test<13, 25, 50, 3>(); end_to_end_test<13, 25, 60, 3>();
	end_to_end_test<13, 60, 50, 3>(); end_to_end_test<13, 60, 60, 3>();
    end_to_end_test<17, 29, 50, 4>(); end_to_end_test<17, 29, 60, 4>(); 
	end_to_end_test<17, 60, 50, 4>(); end_to_end_test<17, 60, 60, 4>();
}
