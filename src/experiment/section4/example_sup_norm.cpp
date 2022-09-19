#include "HEAAN/bootstrap.h"
#include "experiment/rns.h"

#include <iostream>

const double C2 = 5;

template<int LOGN>
void check_sup_norm_bound(){
    std::cout << "Check sup norm is bounded as expected" << std::endl;
    std::cout << "LOGN : " << LOGN << std::endl;

    const int N = 1 << LOGN;

    double max_sup_norm_ratio = 0;
	for(int r = 0; r < 100; ++r) {
		Message<LOGN> z;
		set_random_message(z);
		R_Q<LOGQ, N> pt;
		encode(z, Delta, pt);

		double norm = 0, sup_norm = 0;
		for(int i = 0; i < N; ++i) {
			Z_Q<LOGQ> val = pt[i];
			if(val.is_bigger_than_halfQ())
				val.negate();
			double val_abs_double = (double) val[0];
			norm += val_abs_double * val_abs_double;
			sup_norm = val_abs_double > sup_norm ? val_abs_double : sup_norm;
		}
		norm = sqrt(norm);

        double sup_norm_ratio = sup_norm / (norm / sqrt(N));
        max_sup_norm_ratio = sup_norm_ratio > max_sup_norm_ratio ? sup_norm_ratio : max_sup_norm_ratio;
	}

    double C2_empirical = ceil(max_sup_norm_ratio);
    double C2 = 5;
    std::cout << "max sup_norm / avg_norm ratio : " << max_sup_norm_ratio << std::endl;
    std::cout << "C2(empirical) : " << C2_empirical << std::endl;
    std::cout << "C2 : " << C2 << std::endl;
}

template<int LOGN, int LOGDELTA_boot_tilde, int G>
void measure_sup_norm()
{
    std::cout << "Measuring sup norm on cts input" << std::endl;
    std::cout << "LOGN : " << LOGN << std::endl;
	std::cout << "LOGDELTA_boot_tilde : " << LOGDELTA_boot_tilde << std::endl;
    std::cout << "G : " << G << std::endl;
    const int N = 1 << LOGN;
    const uint64_t Delta_boot_tilde = 1ULL << LOGDELTA_boot_tilde;
    
    // Step 1. Compute constants
    assert((LOGN - 1) % G == 0);
    int D = (LOGN - 1) / G;
    int Diag = pow(3, G);

    double C1 = D * sqrt((H+1) * Diag) / 12.0 * pow(2, 1.0 / (2*D));
    std::cout << "C1 : " << C1 << std::endl;

    // Step 2. Setup CoeffToSlot
    int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z;
    set_random_message(z);

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

    // estimate sup_norm of pt_q / Delta
    double pt_per_Delta_norm = 0;
    double pt_per_Delta_sup_norm = 0;
    for(int i = 0; i < N; ++i) {
        double val = (double) pt_q[i];
        double val_abs_double = std::abs(val) / Delta;

        pt_per_Delta_sup_norm = val_abs_double > pt_per_Delta_sup_norm ? val_abs_double : pt_per_Delta_sup_norm;
        pt_per_Delta_norm += val_abs_double * val_abs_double;
    }
    pt_per_Delta_norm = sqrt(pt_per_Delta_norm);

    std::cout << "sup_norm(pt) bounded : " << pt_per_Delta_norm / sqrt(N) * C2 << std::endl;
    std::cout << "sup_norm(pt) measured : " << pt_per_Delta_sup_norm << std::endl;

    Message<LOGN> z_cts[2], z_cts_exact[2], e[2];
    Message<LOGN+1> e_whole;
    R_Q<LOGQ, N> pt_cts[2];
    R_Q_square<LOGQ,N> ct_cts[2];

    // do cts
    CoeffToSlot<LOGQ,LOGN, LOGDELTA_boot_tilde, G>(ct,s,ct_cts);
    for(int i = 0; i < 2; ++i) {
        HEAAN<LOGQ,N>::dec(ct_cts[i],s,pt_cts[i]);
        decode_log(pt_cts[i],LOGDELTA +(LOGN-1)/G*LOGDELTA_boot_tilde,z_cts[i]);
    }

    // compute desired z_cts
    for(int i = 0; i < N/2; ++i) {
        double val_double = (double) pt[bitReverse(i, LOGN - 1)];
        z_cts_exact[0].r[i] = val_double / Delta;
    }
    for(int i = 0; i < N/2; ++i) {
        double val_double = (double) pt[bitReverse(i, LOGN - 1) + N/2];
        z_cts_exact[1].r[i] = val_double / Delta;
    }

    for(int i = 0; i < 2; ++i)
        sub(z_cts[i], z_cts_exact[i], e[i]);

    aggregate(e[0], e[1], e_whole);

    // estimate sup_norm of e
    double e_norm_expected = C1 / (double) Delta_boot_tilde * pow(N, (1 + 1.0 / (2*D))) * (1 << (LOGq - LOGDELTA));
    double e_norm_measured = norm(e_whole);
    std::cout << "norm(rounding error) expected : " << e_norm_expected << std::endl;
    std::cout << "norm(rounding error) measured : " << e_norm_measured << std::endl;
    std::cout << "sup_norm(rounding error) bounded : " << e_norm_expected / sqrt(N) * C2 << std::endl;
    std::cout << "sup_norm(rounding error) measured : " << sup_norm(e_whole) << std::endl;
}

int main()
{
    //check_sup_norm_bound<9>();
    //check_sup_norm_bound<13>();
    //check_sup_norm_bound<17>();
    measure_sup_norm<9, 50, 2>();
    //measure_sup_norm<13, 50, 3>();
    //measure_sup_norm<17, 50, 4>();
}