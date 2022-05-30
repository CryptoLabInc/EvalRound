#include "HEAAN/bootstrap.h"
#include "experiment/rns.h"

#include <iostream>

template<int LOGN, int G>
void design_delta_tilde()
{
    std::cout << "Designing delta_tilde" << std::endl;
    std::cout << "LOGN : " << LOGN << std::endl;
    std::cout << "G : " << G << std::endl;

    const int N = 1 << LOGN;

    // Step 0. We got C1 as ...
    assert((LOGN - 1) % G == 0);
    int D = (LOGN - 1) / G;
    int Diag = pow(3, G);

    double C1 = D * sqrt((H+1) * Diag) / 12.0 * pow(2, 1.0 / (2*D));
    std::cout << "C1 : " << C1 << std::endl;

    double C2 = 5;
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

int main()
{
    design_delta_tilde<9, 2>();
    design_delta_tilde<17, 4>();
}