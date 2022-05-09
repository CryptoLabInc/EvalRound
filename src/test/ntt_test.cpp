#include <iostream>

#include "HEAAN/advanced/ntt.h"

int main(){
    const int N = 1 << 10;
    NTT<N> ntt(1152921504606877697ULL, 418639631973566421ULL);

    uint64_t a[N];
    for(int i = 0; i < N; ++i) {
        a[i] = i;
    }
    ntt.ntt(a);
    for(int i = 0; i < std::min(N, 10); ++i) {
        std::cout << a[i] << " ";
    }
    std::cout << std::endl;
    ntt.intt(a);
    for(int i = 0; i < std::min(N, 10); ++i) {
        std::cout << a[i] << " ";
    }
}