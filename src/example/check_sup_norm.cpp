#include "HEAAN/bootstrap.h"
#include "experiment/rns_debug.h"

#include <iostream>

int main()
{
    // Step 1. Compute constants
    assert((LOGN - 1) % G == 0);
    int D = (LOGN - 1) / G;
    int Diag = pow(3, G);

    double C1 = D * sqrt((H+1) * Diag) / 12.0 * pow(2, 1.0 / (2*D));
    std::cout << "C1 : " << C1 << std::endl;

    // Step 2. Do CoeffToSlot
    int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z;
    set_random_message(z);

    // get z_amb
    R_Q<LOGq, N> pt_q;
    R_Q_square<LOGq,N> ct_q;
    R_Q_square<LOGQ,N> ct;
    R_Q<LOGQ, N> pt;
    Message<LOGN> z_amb;
    encode(z,Delta,pt_q);
    HEAAN<LOGq,N>::enc(pt_q,s,ct_q);
    resize(ct_q, ct); // size up to Rq -> RQ
    HEAAN<LOGQ,N>::dec(ct,s, pt);
    decode(pt, Delta, z_amb);

    Message<LOGN> z_cts[2], z_cts_exact[2], e[2];
    Message<LOGN+1> e_whole;
    R_Q<LOGQ, N> pt_cts[2];
    R_Q_square<LOGQ,N> ct_cts[2];

    // do cts, gather splitted z_cts
    CoeffToSlot<LOGQ,LOGN, LOGDELTA_boot_tilde, G>(ct,s,ct_cts);
    for(int i = 0; i < 2; ++i) {
        HEAAN<LOGQ,N>::dec(ct_cts[i],s,pt_cts[i]);
        decode_log(pt_cts[i],LOGDELTA +(LOGN-1)/G*LOGDELTA_boot_tilde,z_cts[i]);
    }

    // compute desired z_cts
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

        z_cts_exact[0].r[i] = valr_double;
        z_cts_exact[1].r[i] = vali_double;
    }
    for(int i = 0; i < 2; ++i)
        sub(z_cts[i], z_cts_exact[i], e[i]);
    for(int i = 0; i < N; ++i) {
        e_whole.r[i] = i < N/2 ? e[0].r[i] : e[1].r[i];
        e_whole.i[i] = 0;
    }

    // estimate sup_norm of pt_q / Delta
    double pt_per_Delta_norm = 0;
    double pt_per_Delta_sup_norm = 0;
    for(int i = 0; i < N; ++i) {
        Z_Q<LOGq> val = pt_q[i];
        if(val.is_bigger_than_halfQ())
            val.negate();
        double val_abs_double = (double) val[0] / Delta;

        pt_per_Delta_sup_norm = val_abs_double > pt_per_Delta_sup_norm ? val_abs_double : pt_per_Delta_sup_norm;
        pt_per_Delta_norm += val_abs_double * val_abs_double;
    }
    pt_per_Delta_norm = sqrt(pt_per_Delta_norm);

    std::cout << "sup_norm(pt) bounded : " << pt_per_Delta_norm / sqrt(N) * 5 << std::endl;
    std::cout << "sup_norm(pt) measured : " << pt_per_Delta_sup_norm << std::endl;

    // estimate sup_norm of e
    double e_norm_expected = C1 / (double) Delta_boot_tilde * pow(N, (1 + 1.0 / (2*D))) * (1 << (LOGq - LOGDELTA));
    double e_norm_measured = norm(e_whole);
    std::cout << "norm(rounding error) expected : " << e_norm_expected << std::endl;
    std::cout << "norm(rounding error) measured : " << e_norm_measured << std::endl;
    std::cout << "sup_norm(rounding error) bounded : " << e_norm_expected / sqrt(N) * 5 << std::endl;
    std::cout << "sup_norm(rounding error) measured : " << sup_norm(e_whole) << std::endl;

/*    // evaluate eval_mod on z_cts_res and z_cts_exact
    
    const int LOGQ_after_cts = LOGQ-(LOGN-1)/G*LOGDELTA_boot_tilde;
	R_Q_square<LOGQ_after_cts,N> ct_ctsrs[2];
	for(int i=0;i<2;i++)
	    RS<LOGQ,LOGQ_after_cts,N>(ct_cts[i],ct_ctsrs[i]);

    R_Q<LOGQ_after_cts, N> pt_cts_exact;
    R_Q_square<LOGQ_after_cts, N> ct_cts_exact;

  	const int LOGQ_after_evalmod = LOGQ_after_cts - 12*LOGDELTA_boot;
    R_Q_square<LOGQ_after_evalmod, N> ct_evalmod, ct_evalmod_exact;
    R_Q<LOGQ_after_evalmod, N> pt_evalmod, pt_evalmod_exact;
    Message<LOGN> z_evalmod, z_evalmod_exact;

    for(int i = 0; i < N/2; ++i) {
        z_cts_exact.i[i] = 0;
    }
    encode(z_cts_exact,Delta,pt_cts_exact);
    HEAAN<LOGQ_after_cts,N>::enc(pt_cts_exact,s,ct_cts_exact);
    EvalMod<LOGQ_after_cts,N,LOGDELTA_boot,K>(ct_ctsrs[0],s,ct_evalmod);
    EvalMod<LOGQ_after_cts,N,LOGDELTA_boot,K>(ct_cts_exact,s,ct_evalmod_exact);
    HEAAN<LOGQ_after_evalmod,N>::dec(ct_evalmod,s,pt_evalmod);
    HEAAN<LOGQ_after_evalmod,N>::dec(ct_evalmod_exact,s,pt_evalmod_exact);
	decode_log(pt_evalmod, LOGDELTA, z_evalmod);
    decode_log(pt_evalmod_exact, LOGDELTA, z_evalmod_exact);
    sub(z_evalmod, z_evalmod_exact, e);

    // estimate sup_norm of final error
    std::cout << "sup_norm(result error) bounded : " << (pt_per_Delta_norm / sqrt(N) * 5 + e_norm_expected / sqrt(N/2) * 5 ) << std::endl;
    std::cout << "sup_norm(result error) measured : " << sup_norm(e) << std::endl;*/
}