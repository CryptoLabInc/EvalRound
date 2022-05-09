#include <cstdint>

template <int N>
class NTT {
    uint64_t q;
    uint64_t psi_rev[N];
    uint64_t inv_psi_rev[N];
    uint64_t inv_N;

    NTT(uint64_t q, uint64_t psi);

    void ntt(uint64_t a[N]);
    void intt(uint64_t a[N]);

    private:
        uint64_t power_mod(uint64_t a, int n);
        uint64_t inv_mod(uint64_t a);
};

inline int log(int N) {
    int i = 0, n = 1;
    while(n < N) {
        i++;
        n *= 2;
    }
    return i;
}

inline uint32_t bitReverse32(uint32_t x) {
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return ((x >> 16) | (x << 16));
}

inline uint32_t bitReverseN(uint32_t x, int LOGN) {
    return bitReverse32(x) >> (LOGN - 1);
}

template<int N>
NTT<N>::NTT(uint64_t q, uint64_t psi) {
    this->q = q;
    uint64_t psi[N];
    psi[0] = 1;
    for(int i = 1; i < N; ++i) {
        psi[i] = psi
    }
    for(int i = 0; i < N; ++i) {
        this->psi_rev[i] = 
    }
    this->psi_rev
}

template<int N>
uint64_t NTT<N>::power_mod(uint64_t a, uint64_t n) {
    assert(pow >= 0);
    Z_P<P> temp(1);
    for(int i = 63; i >= 0; --i){
        int digit = (pow >> i) & 1;
        temp = temp * temp;
        if(digit == 1)
            temp = temp * (*this);
    }
    return temp;
}

template<int N>
uint64_t NTT<N>::inv_mod(uint64_t q, uint64_t psi);