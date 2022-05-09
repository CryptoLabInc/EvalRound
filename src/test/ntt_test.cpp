#include <iostream>

#include "HEAAN/advanced/ntt.h"

void ntt_test() {
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

void mod_test() {
    uint64_t a = (uint64_t) 1 << 63, b = ((uint64_t) 1 << 63) + 1, q = 100;
    std::cout << add_mod(a, b, q) << std::endl;
    std::cout << sub_mod(a, b, q) << std::endl;
    std::cout << mul_mod(a, b, q) << std::endl;
    std::cout << power_mod(a, 5, q) << std::endl;
}

int main(){
    ntt_test();
}