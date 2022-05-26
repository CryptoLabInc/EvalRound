#include "experiment/lattigo.h"
#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_bootstrap.h"

#include <iostream>
#include <cmath>

int main()
{
    double q = 1ULL << LOGq;
    double z_amb_norm_expected = sqrt(N/2) / Delta * sqrt((H+1.0)/12*q*q*N);
    int G[4] = {4, 4, 4, 3};
    int Diag[4];
    for(int i = 0; i < 4; ++i)
        Diag[i] = pow(3, G[i]);
    double p_A[4];
    for(int i = 0; i < 4; ++i)
        p_A[i] = sqrt(Diag[i]*N / 12.0) / Delta_boot_tilde;
    double U0_norm_one = sqrt(2) / pow(N, 1.0 / (LOGN - 1));
    double A_norm[4];
    for(int i = 0; i < 4; ++i)
        A_norm[i] = pow(U0_norm_one, G[i]);
    double theorem1_const = 0;
    for(int i = 0; i < 4; ++i) {
        double norm_const = 1;
        for(int j = 0; j < 4; ++j) {
            if(i == j)
                continue;
            norm_const *= A_norm[j];
        }
        theorem1_const += p_A[i] * norm_const;
    }
    double e_norm_expected = z_amb_norm_expected * theorem1_const * 2;

    std::cout << e_norm_expected * 5 / sqrt(N) << std::endl;
    double eps = pow(2.0, -15);
    std::cout << eps * (1 << (LOGq - LOGDELTA)) / 2 << std::endl;
}
