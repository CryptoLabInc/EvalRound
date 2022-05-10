#include <iostream>

#include "HEAAN/advanced/crt.h"

int main(){
    const int L = 5;
    const int LOGQ = 64*5;
    uint64_t q[L]={1152921504606847009ULL,1152921504606847067ULL,1152921504606847081ULL,
        1152921504606847123ULL,1152921504606847127ULL};

    CRT<L, LOGQ> crt(q);

    uint64_t a_data[L] = {2496847982066394524ULL,6179808641469508490ULL,4613858499374131512ULL,
14039320798214412288ULL,10083004203411ULL};

    Z_Q<LOGQ> a;
    for(int i = 0; i < L; ++i){
        a[i]= a_data[i];
    }

    uint64_t a_rns[L] = {753612675309782689ULL,8006693014835917ULL,16048976969132353ULL,376433015983029137ULL,
148789351586435012ULL};

    uint64_t a_rns_measured[L] = {0};
    crt.icrt(a, a_rns_measured);
    for(int i = 0; i < L; ++i)
        std::cout << a_rns[i] - a_rns_measured[i] << " ";
    std::cout << std::endl;

    for(int i = 0; i < L; ++i)
        std::cout << a_rns[i] << " ";
    std::cout << std::endl;

    /*Z_Q<LOGQ> a_measured;
    crt.crt(a_rns, a_measured);
    for(int i = 0; i < L; ++i)
        std::cout << a[i] - a_measured[i] << " ";
    std::cout << std::endl;*/
}