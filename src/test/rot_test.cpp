#include "experiment/big.h"

#include "util/util.h"

#include <iostream>

int main(){
    int s[N], s_rot[N] = {0}, s_rot_new[N] = {0};
    HEAAN<LOGQ,N>::keygen(4,s);
    for(int r = 0; r < 10; ++r) {
        std::cout << std::endl;
        int rotation = rand() % (N/2);
        Timer t1("rot");
        rot<N>(s, s_rot, rotation);
        t1.stop();
        Timer t2("rot_new");
        rot_new<N>(s, s_rot_new, rotation);
        t2.stop();
        for(int i = 0; i < N; ++i) {
            if(s_rot[i] != s_rot_new[i])
                std::cout << i << " " << s_rot[i] << " " << s_rot_new[i] << std::endl;
        }
    }
}