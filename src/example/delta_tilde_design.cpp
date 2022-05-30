#include "HEAAN/bootstrap.h"
#include "experiment/rns_debug.h"

#include <iostream>

int main()
{
    // Step 0. We got C1 as ...
    assert((LOGN - 1) % G == 0);
    int D = (LOGN - 1) / G;
    int Diag = pow(3, G);

    double C1 = D * sqrt((H+1) * Diag) / 12.0 * pow(2, 1.0 / (2*D));
    std::cout << "C1 : " << C1 << std::endl;

    // Step 1. Empirically get C2
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

		double avg_norm = norm * (1 / sqrt(N));
		double ratio = sup_norm / avg_norm;
        max_sup_norm_ratio = ratio > max_sup_norm_ratio ? ratio : max_sup_norm_ratio;
	}

    double C2_empirical = ceil(max_sup_norm_ratio);
    double C2 = 5;
    std::cout << "max sup_norm / avg_norm ratio : " << max_sup_norm_ratio << std::endl;
    std::cout << "C2(empirical) : " << C2_empirical << std::endl;
    std::cout << "C2 : " << C2 << std::endl;

    // Step 2. Compute Delta given value \simeq rounding error
    double Delta_tilde_sim = C1 * pow(N, (1 + 1.0 / (2*D))) * (1 << (LOGq - LOGDELTA));
    int LOGDELTA_tilde_sim = ceil(log2(Delta_tilde_sim));
    std::cout << "Delta_tilde (value is similar to rounding error) : " << LOGDELTA_tilde_sim << std::endl;


    // Step 3. Compute Delta when value + rounding error < epsilon
    double epsilon = 1;
    std::cout << "given epsilon (max distance from multiple of q / Delta) : " << epsilon << std::endl;
    double tight_ratio = 1.0 / (epsilon * sqrt(N) / C2 - 1);
    int Delta_tilde_tight = Delta_tilde_sim * tight_ratio;
    int LOGDELTA_tilde_tight = ceil(log2(Delta_tilde_tight));
    std::cout << "Delta_tilde (value + rounding error < epsilon = 1) : " << LOGDELTA_tilde_tight << std::endl;
}