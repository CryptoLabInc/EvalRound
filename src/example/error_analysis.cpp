#include "experiment/rns_debug.h"
#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"

#include <iostream>
#include <cmath>

int main()
{
    // Step 1. Compute constants
    assert((LOGN - 1) % G == 0);
    int D = (LOGN - 1) / G;
    int Diag = pow(3, G);

    double C1 = D * sqrt((H+1) * Diag) / 24 * pow(2, 1.0 / (2*D)) * sqrt(2);
    std::cout << "C1 : " << C1 << std::endl;

    // Step 2. Measure norm(z_amb)
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

    // compare
    double q = 1ULL << LOGq;
    double z_amb_norm_expected = sqrt(N/2) / Delta * sqrt((H+1)/12*q*q*N);;
    double z_amb_norm_measured = norm(z_amb);
    std::cout << "norm(z_amb) expected : " << z_amb_norm_expected << std::endl;
    std::cout << "norm(z_amb) measured : " << z_amb_norm_measured << std::endl;


    // Step 3. Compare norm(rounding error of cts)
    Message<LOGN> z_cts[2], z_cts_res, z_cts_exact, e;
    R_Q<LOGQ, N> pt_cts[2];
    R_Q_square<LOGQ,N> ct_cts[2];

    // do cts, gather splitted z_cts
    CoeffToSlot<LOGQ,LOGN, LOGDELTA_boot_tilde, G>(ct,s,ct_cts);
    for(int i = 0; i < 2; ++i) {
        HEAAN<LOGQ,N>::dec(ct_cts[i],s,pt_cts[i]);
        decode_log(pt_cts[i],LOGDELTA +(LOGN-1)/G*LOGDELTA_boot_tilde,z_cts[i]);
    }
    for(int i = 0; i < N/2; ++i) {
        z_cts_res.r[i] = z_cts[0].r[i];
        z_cts_res.i[i] = z_cts[1].r[i];
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

        z_cts_exact.r[i] = valr_double;
        z_cts_exact.i[i] = vali_double;
    }

    // compare
    sub(z_cts_res, z_cts_exact, e);
    double e_norm_expected = C1 / (double) Delta_boot_tilde * pow(N, (1 + 1.0 / (2*D))) * (1 << (LOGq - LOGDELTA));
    double e_norm_measured = norm(e);
    std::cout << "norm(e) expected : " << e_norm_expected << std::endl;
    std::cout << "norm(e) measured : " << e_norm_measured << std::endl;

    // sanity check
    {
        double q = 1ULL << LOGq;
        double z_amb_norm_expected = sqrt(N/2) / Delta * sqrt((H+1.0)/12*q*q*N);
        double p_U0 = sqrt(Diag*N / 12.0) / Delta_boot_tilde;
        double U0_norm = pow(sqrt(2) / pow(N, 1.0 / (LOGN - 1)), G);
        double e_norm_expected = z_amb_norm_expected * p_U0 * D * pow(U0_norm, D - 1) * sqrt(2);
        std::cout << "norm(e) expected (sanity check) : " << e_norm_expected << std::endl;
    }
}