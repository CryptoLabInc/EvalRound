#include "experiment/big.h"

#include "util/util.h"

#include <iostream>

int main(){
    int s[N], s_rot_old[N] = {0}, s_rot[N] = {0};
    HEAAN<LOGQ,N>::keygen(4,s);
    for(int r = 0; r < 10; ++r) {
        std::cout << std::endl;
        int rotation = rand() % (N/2);
        Timer t1("rot_old");
        rot_old<N>(s, s_rot_old, rotation);
        t1.stop();
        Timer t2("rot");
        rot<N>(s, s_rot, rotation);
        t2.stop();
        for(int i = 0; i < N; ++i) {
            if(s_rot_old[i] != s_rot[i])
                std::cout << i << " " << s_rot_old[i] << " " << s_rot[i] << std::endl;
        }
    }
}