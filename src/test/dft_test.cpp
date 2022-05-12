#include "experiment/simple.h"

#include "HEAAN/DFT.h"

void ifft_test(){
    Message<LOGN> z;
    double m_idft[N], m_ifft[N];
    set_test_rounded_message(z, Delta);
    idft<N>(z.r, z.i, m_idft);
    ifft<LOGN>(z.r, z.i, m_ifft);
    print_array<N>(m_idft);
    print_array<N>(m_ifft);
}

void run_test(){
    Message<LOGN> z, z_res;
    double m[N];

    set_random_message(z);
    ifft<LOGN>(z.r, z.i, m);
    fft<LOGN>(m, z_res.r, z_res.i);
    print("zr", z);
    print("zr_res", z_res);
}


int main()
{
  ifft_test();
  run_test();
}
