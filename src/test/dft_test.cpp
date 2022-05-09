#include "experiment/simple.h"

#include "HEAAN/DFT.h"

void ifft_test(){
    Message<LOGN> z;
    double m_idft[N], m_ifft[N];
    set_test_rounded_message(z, Delta);
    idft<N>(z.r, z.i, m_idft);
    ifft<LOGN>(z.r, z.i, m_ifft);
    for(int i = 0; i < std::min(N, 10); ++i) {
        std::cout << i << " " << m_idft[i] << " " << m_ifft[i] << std::endl;
    }
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
