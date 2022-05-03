#include "SIMPLE/test.h"

#include "HEAAN/DFT.h"

void fft_test(){
    double m[N];
    double zr_dft[N/2], zi_dft[N/2];
    double zr_fft[N/2], zi_fft[N/2];
    double er[N/2], ei[N/2];

    set_random_message(m, m+N/2);
    dft<N>(m, zr_dft, zi_dft);
    fft<LOGN>(m, zr_fft, zi_fft);
    print("zr_dft", zr_dft);
    print("zr_fft", zr_fft);
}

void fft_run_test(){
    double m[N];
    double zr[N/2], zi[N/2];
    double zr_res[N/2], zi_res[N/2];

    set_random_message(zr, zi);
    ifft<LOGN>(zr, zi, m);
    fft<LOGN>(m, zr_res, zi_res);
    print("zr", zr);
    print("zr_res", zr_res);
}

void idft_test(){
    double zr[N/2], zi[N/2];
    double m_idft[N], m_ifft[N];
    set_test_rounded_message(zr, zi);
    idft<N>(zr, zi, m_idft);
    ifft<LOGN>(zr, zi, m_ifft);
    print_pt("m_idft", m_idft);
    print_pt("m_ifft", m_ifft);
}
int main()
{
  fft_test();
  idft_test();
}
