#include "HEAAN/bootstrap.h"
#include "experiment/rns_debug.h"

#include <iostream>

void measure_b_evalmod()
{
    int s[N], s_sq[N];
    HEAAN<LOGQ,N>::keygen(H,s);

    Message<LOGN> z, z_res;

    int eps_num = N / 2 / 49 / 2;
    double q_per_Delta = pow(2, 10);
	for(int i = 0; i < 49*2*eps_num; ++i) {
        double val = (i%(2*eps_num) - eps_num) / (double) eps_num;
        int I = q_per_Delta * ((i / (2*eps_num)) - 24);
        z.r[i] = val + I;
        z.i[i] = 0;
        z_res.r[i] = val;
        z_res.i[i] = 0;
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
    sub(z_res, z_evalmod, e);
    print("z", z);
    print("z_res", z_res);
    print("z_evalmod", z_evalmod);
    print("e", e);
    
    Message<LOGN> er, ez;
    for(int i = 0; i < N/2; ++i) {
        er.r[i] = e.r[i];
        ez.i[i] = e.i[i];
    }
    std::cout << "sup_norm : " << sup_norm(e) << std::endl;
    std::cout << "sup_norm : " << sup_norm(er) << " " << sup_norm(ez) << std::endl;
}

int main()
{
	measure_b_evalmod();
}