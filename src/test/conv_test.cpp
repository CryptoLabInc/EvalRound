#include "experiment/test.h"
#include "HEAAN/arith/conv.h"

#include <iostream>

int main(){    
    Message<LOGN> a, b, c_orig, c;
    R_Q<LOGQ, N> pt_a, pt_b, pt_c_orig, pt_c;

    set_random_message(a);
    set_random_message(b);
    encode(a, Delta, pt_a);
    encode(b, Delta, pt_b);
    
    mul(a, b, c_orig);
    print("c_orig", c_orig);

    conv(pt_a, pt_b, pt_c);
    decode_log(pt_c, LOGDELTA*2, c);
    print("c", c);
}