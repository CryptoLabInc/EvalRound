#include "impl/message.h"
#include "impl/simple_plaintext.h"

#include "HEAAN/R_Q.h"
#include "HEAAN/DFT.h"

#include <random>

namespace {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1.0, 1.0);
}

template <class Iter>
void sampleUniform(Iter begin, Iter end) {
    for(Iter it = begin; it != end; ++it) {
        *it = dist(gen);
    }
}

template <int LOGN>
void set_test_message(Message<LOGN> &z) {
    const int N = 1 << LOGN;
    for(int i = 0; i < N/2; ++i) {
        double x = (double) (i) / (N/2) * PI;
        z.r[i] = cos(x) / sqrt(2);
        z.i[i] = sin(x) / sqrt(2);
    }
}

template <int LOGN>
void set_random_message(Message<LOGN> &z) {
    const int N = 1 << LOGN;
    sampleUniform(z.r, z.r+N/2);
    sampleUniform(z.i, z.i+N/2);
}

template <int LOGN>
void set_test_rounded_message(Message<LOGN> &z, uint64_t Delta) {
    set_test_message(z);
    SimplePlaintext<LOGN> pt;
    encode(z, Delta, pt);
    decode(pt, Delta, z);
}

template <int LOGN, int K>
void set_test_matrix(Message<LOGN> A[K]) {
    const int N = 1 << LOGN;
    for(int k = 0; k < K; ++k) {
        for(int i = 0; i < N / 2; ++i) {
            double x = (double) (k*i) / (N/2) * PI;
            A[k].r[i] = cos(x) / sqrt(2);
            A[k].r[i] = sin(x) / sqrt(2);
        }
    }
}

template <int LOGN, int K>
void set_random_matrix(Message<LOGN> A[K]) {
    for(int k = 0; k < K; ++k) {
        set_random_message(A[k]);
    }
}

template <int LOGN>
void set_test_U0_matrix(SparseDiagonal<(1<<(LOGN-1)),3> U0r[LOGN-1],
	            SparseDiagonal<(1<<(LOGN-1)),3> U0i[LOGN-1]) {
    splitU0NR<LOGN>(U0r, U0i);
    
    for (int n = 0; n < LOGN - 1; n++) {
		U0r[n].transpose();
		U0i[n].transpose();
		U0i[n].negate();
	}
}

template <int LOGQ, int N>
void set_test_plaintext(R_Q<LOGQ, N> &A){
    for(int i = 0; i < N; ++i) {
        A[i].setzero();
        A[i][0] = i;
    }
}