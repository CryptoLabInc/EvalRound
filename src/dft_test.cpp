#include "test.h"

#include "HEAAN/DFT.h"

void fft_test(){
    double m[N];
    double zr_dft[N/2], zi_dft[N/2];
    double zr_fft[N/2], zi_fft[N/2];
    double er[N/2], ei[N/2];

    set_test_message(m);
    dft<N>(m, zr_dft, zi_dft);
    fft<N>(m, zr_fft, zi_fft);
    sub(zr_dft, zr_fft, er);
    sub(zi_dft, zi_fft, ei);
    print("er", er);
    print("ei", ei);
}

void fft_run_test(){
    double m[N];
    double zr_fft[N/2], zi_fft[N/2];

    set_test_message(m);
    fft<N>(m, zr_fft, zi_fft);
    print("zr", zr_fft);
}

void idft_correct_test(){
    double m[N], m_res[N];
    double zr[N/2], zi[N/2];
    
    set_test_message(m);
    dft_orig<N>(m, zr, zi);
    idft_orig<N>(zr, zi, m_res);
    print_pt("m", m);
    print_pt("m_res", m_res);
}

void idft_test(){
    double zr[N/2], zi[N/2];
    double m_idft[N], m_ifft[N];
    
    set_test_message(zr, zi);
    print("zr", zr);
    print("zi", zi);
    idft_orig<N>(zr, zi, m_idft);
    ifft_orig<N>(zr, zi, m_ifft);
    print_pt("m_idft", m_idft);
    print_pt("m_ifft", m_ifft);

    idft<N>(zr, zi, m_idft);
    ifft<N>(zr, zi, m_ifft);
    print_pt("m_idft", m_idft);
    print_pt("m_ifft", m_ifft);
}
int main()
{
    idft_correct_test();
}