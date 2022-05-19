#include "impl/message.h"
#include "impl/simple_plaintext.h"

#include "HEAAN/R_Q.h"
#include "HEAAN/DFT.h"

#include <cmath>
#include <random>

namespace {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    double range = 1. / sqrt(2);
    std::uniform_real_distribution<> dist(-range, range);
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
        z.r[i] = cos(x);
        z.i[i] = sin(x);
    }
}

template <int LOGN>
void print_max_error(const Message<LOGN> &z, const Message<LOGN> &z_out) {
    const int N = 1 << LOGN;
	double max_err_real = 0.0, max_err_imag = 0.0;
	for(int i = 0; i < N/2; ++i) {
		max_err_real = std::max(max_err_real, std::fabs(z.r[i]-z_out.r[i]));
		max_err_imag = std::max(max_err_imag, std::fabs(z.i[i]-z_out.i[i]));
	}
	std::cout << "LOG2 MAX ERROR (REAL, IMAG) = [" << std::log2(max_err_real) << ", " << std::log2(max_err_imag) << "]" << std::endl;
}

template <int LOGN>
void set_random_message(Message<LOGN> &z) {
    const int N = 1 << LOGN;
    sampleUniform(z.r, z.r+N/2);
    sampleUniform(z.i, z.i+N/2);
}

template <int LOGN>
void set_evalmod_message(Message<LOGN> &z) {
    const int N = 1 << LOGN;
	for(int i = 0; i < N/2; i++) {
	    z.r[i] = 0.0001 + i % 24;
		z.i[i] = 0;
	}
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
