#include <iostream>

#include "HEAAN/advanced/ntt.h"

int main(){
    const int N = 1 << 10;
    NTT<N> ntt(1152921504606877697ULL);

    uint64_t a[N];
    for(int i = 0; i < N; ++i) {
        a[i] = i;
    }
    ntt.ntt(a);
    print_array<N>(a);
    ntt.intt(a);
    print_array<N>(a);
}