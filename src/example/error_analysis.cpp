#include "experiment/rns_debug.h"

#include <iostream>

int main()
{
    int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z;
    set_random_message(z);

    R_Q<LOGq, N> pt_q;
    R_Q_square<LOGq,N> ct_q;
    R_Q_square<LOGQ,N> ct_Q;
    R_Q<LOGQ, N> pt_Q;
    Message<LOGN> z_amb;
    encode(z,Delta,pt_q);
    HEAAN<LOGq,N>::enc(pt_q,s,ct_q);
    resize(ct_q, ct_Q); // size up to Rq -> RQ
    HEAAN<LOGQ,N>::dec(ct_Q,s, pt_Q);
    decode(pt_Q, Delta, z_amb);

    double q = 1ULL << LOGq;
    double norm_expected = sqrt(N/2) / Delta * sqrt((H+1)/12*q*q*N);
    double norm_measured = norm(z_amb);
    std::cout << "norm(z_amb) expected : " << norm_expected << std::endl;
    std::cout << "norm(z_amb) measured : " << norm_measured << std::endl;


}