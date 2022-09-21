#include "HEAAN/bootstrap.h"
#include "experiment/rns.h"

#include <iostream>

template<int LOGN>
void measure_b_evalmod()
{
    std::cout << "Measuring error bound of evalmod result" << std::endl;
    std::cout << "LOGN : " << LOGN << std::endl;

    const int N = 1 << LOGN;

    int s[N], s_sq[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z, z_evalmod_exact;

    int eps_num = (N / 2 / 49 - 1) / 2;
    double q_per_Delta = pow(2, 10);
    int idx = 0;
    for(int i = -24; i <= 24; ++i) {
        for(int j = -eps_num; j <= eps_num; ++j) {
            z.r[idx] = i * q_per_Delta + j / (double) eps_num;
            z_evalmod_exact.r[idx] = j / (double) eps_num;
            idx++;
        }
    }

    R_Q<LOGQ, N> pt;
    R_Q_square<LOGQ, N> ct;
    encode(z, Delta, pt);
    HEAAN<LOGQ,N>::enc(pt,s,ct);

	// EvalMod
	const int LOGQ_after_evalmod = LOGQ - 12*LOGDELTA_boot;
	R_Q_square<LOGQ_after_evalmod,N> ct_evalmod;
	for(int i=0;i<2;i++)
	    EvalMod<LOGQ,N,LOGDELTA_boot,K>(ct,s,ct_evalmod);
    
	R_Q<LOGQ_after_evalmod,N> pt_evalmod;
	HEAAN<LOGQ_after_evalmod,N>::dec(ct_evalmod,s,pt_evalmod);
    Message<LOGN> z_evalmod;
	decode_log(pt_evalmod,LOGDELTA,z_evalmod);

    Message<LOGN> e;
    sub(z_evalmod_exact, z_evalmod, e);
    
    for(int i = 0; i < N/2; ++i) {
        e.i[i] = 0;
    }
    std::cout << "sup_norm : " << sup_norm(e) << std::endl;
}

int main()
{
	measure_b_evalmod<9>();
    measure_b_evalmod<13>();
    measure_b_evalmod<17>();
}