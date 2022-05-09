#include <iostream>

#include "HEAAN/advanced/ntt.h"

int main(){
    const int N = 1 << 10;
    NTT<N>(1152921504606877697ULL, 418639631973566421ULL);
}