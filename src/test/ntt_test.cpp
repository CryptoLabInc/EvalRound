#include "HEAAN/arith/conv/ntt.h"

#include <iostream>

int main(){
    const int N = 1 << 10;
    NTT<N> ntt(1152921504606877697ULL);

    uint64_t a[N];
    for(int i = 0; i < N; ++i) {
        a[i] = i;
    }
    ntt.ntt(a);
    ntt.intt(a);
    print_array<N>(a);
}