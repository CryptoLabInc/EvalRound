#pragma once

#include "HEAAN/matrix.h"

#include <iostream>
#include <string>
#include <cmath>

template <int LOGN>
struct Message {
    double r[1 << (LOGN - 1)] = {0};
    double i[1 << (LOGN - 1)] = {0};
};

template <int LOGN, int K>
void matrix_vector_product(
    const Message<LOGN> &z,
    const Message<LOGN> A[K],
    Message<LOGN> &Az) {
    const int N = 1 << LOGN;
    for(int i = 0; i < N/2; ++i) {
        double sumr = 0, sumi = 0;
        for(int k = 0; k < K; ++k) {
            int j = (i+k) % (N/2);
            sumr += A[k].r[i] * z.r[j] - A[k].i[i] * z.i[j];
            sumi += A[k].r[i] * z.i[j] + A[k].i[i] * z.r[j];
        }
        Az.r[i] = sumr;
        Az.i[i] = sumi;
    }
}

template <int LOGN>
void matrix_vector_product(
    const Message<LOGN> &z,
    SparseDiagonal<(1<<(LOGN-1)),3> &Ar,
	SparseDiagonal<(1<<(LOGN-1)),3> &Ai,
    Message<LOGN> &Az) {
    const int N = 1 << LOGN;
    for(int i = 0; i < N/2; ++i) {
        double sumr = 0, sumi = 0;
        for(int k = 0; k < 3; ++k) {
            int jr = (i+Ar.off[k]) % (N/2);
            int ji = (i+Ai.off[k]) % (N/2);
            sumr += Ar.vec[k][i] * z.r[jr] - Ai.vec[k][i] * z.i[ji];
            sumi += Ar.vec[k][i] * z.i[ji] + Ai.vec[k][i] * z.r[jr];
        }
        Az.r[i] = sumr;
        Az.i[i] = sumi;
    }
}

template <int LOGN>
void add(const Message<LOGN> &z1, const Message<LOGN> &z2, Message<LOGN> &z3){
	for(int i=0; i< (1 << (LOGN - 1)); i++){
		z3.r[i] = z1.r[i] + z2.r[i];
   		z3.i[i] = z1.i[i] + z2.i[i];
	}
}

template <int LOGN>
void sub(const Message<LOGN> &z1, const Message<LOGN> &z2, Message<LOGN> &z3){
	for(int i=0; i< (1 << (LOGN - 1)); i++){
   		z3.r[i] = z1.r[i] - z2.r[i];
        z3.i[i] = z1.i[i] - z2.i[i];
	}
}

template <int LOGN>
double square_sum(const Message<LOGN> &z){
    double sum = 0;
    for(int i = 0; i < (1 << (LOGN - 1)); ++i) {
        sum += z.r[i] * z.r[i] + z.i[i] * z.i[i];
    }
    return sum;
}

template <int LOGN>
double norm(const Message<LOGN> &z){
    return sqrt(square_sum(z));
}

template <int LOGN>
void print(const std::string name, const Message<LOGN> &z){
	std::cout << "Message " << name << std::endl;
	for(int i=0; i <std::min((1 << (LOGN - 1)), 10); ++i) {
		std::cout << "(" << z.r[i] << ", " << z.i[i] << ") ";
	}
	std::cout << std::endl;
}