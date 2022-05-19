#include <iostream>

#include "HEAAN/advanced/crt.h"

int main(){
    const int L = 5;
    const int LOGQ = 64*5;
    uint64_t q[L]={1152921504606847009ULL,1152921504606847067ULL,1152921504606847081ULL,
        1152921504606847123ULL,1152921504606847127ULL};

    CRT<L, LOGQ> crt(q);

    uint64_t a_data[5] = {2496847982066394524ULL,6179808641469508490ULL,4613858499374131512ULL,
14039320798214412288ULL, 0ULL};

    Z_Q<LOGQ> a;
    for(int i = 0; i < 5; ++i){
        a[i]= a_data[i];
    }

    uint64_t a_rns[L] = {171973778066475490ULL,
    441363031084368163ULL,
    885382643671141350ULL,
    1019151257649711462ULL,
    1052569354770772635ULL};

    uint64_t a_rns_measured[L] = {0};
    crt.icrt(a, a_rns_measured);
    for(int i = 0; i < L; ++i)
        std::cout << a_rns[i] << " " << a_rns_measured[i] << " ";
    std::cout << std::endl;

    Z_Q<LOGQ> a_measured;
    crt.crt(a_rns, a_measured);
    for(int i = 0; i < 5; ++i) {
        std::cout << a[i] - a_measured[i] << " ";
    }
    std::cout << std::endl;
}