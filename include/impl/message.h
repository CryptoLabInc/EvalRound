#pragma once

#include "HEAAN/matrix.h"
#include "HEAAN/DFT.h"
#include "util/util.h"

#include <iostream>
#include <string>
#include <cmath>

template <int LOGN>
struct Message {
    double r[1 << (LOGN - 1)];
    double i[1 << (LOGN - 1)];

    Message<LOGN>() {
        for(int i = 0; i < 1 << (LOGN - 1); ++i){
            this->r[i] = 0;
            this->i[i] = 0;
        }
    }

    Message<LOGN>(const Message<LOGN> &other) {
        for(int i = 0; i < 1 << (LOGN - 1); ++i){
            this->r[i] = other.r[i];
            this->i[i] = other.i[i];
        }
    }
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

template <int LOGN, int S>
void matrix_vector_product(
    const Message<LOGN> &z,
    SparseDiagonal<(1<<(LOGN-1)),S> &Ar,
	SparseDiagonal<(1<<(LOGN-1)),S> &Ai,
    Message<LOGN> &Az) {
    const int N = 1 << LOGN;
    for(int i = 0; i < N/2; ++i) {
        double sumr = 0, sumi = 0;
        for(int k = 0; k < S; ++k) {
            if(!Ar.zero[k] || !Ai.zero[k]){
                int jr = (i+Ar.off[k]) % (N/2);
                int ji = (i+Ai.off[k]) % (N/2);
                sumr += Ar.vec[k][i] * z.r[jr] - Ai.vec[k][i] * z.i[ji];
                sumi += Ar.vec[k][i] * z.i[ji] + Ai.vec[k][i] * z.r[jr];
            }
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
void mul(const Message<LOGN> &z1, const Message<LOGN> &z2, Message<LOGN> &z3){
	for(int i=0; i< (1 << (LOGN - 1)); i++){
   		z3.r[i] = z1.r[i] * z2.r[i] - z1.i[i] * z2.i[i];
        z3.i[i] = z1.r[i] * z2.i[i] + z1.i[i] * z2.r[i];
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

template<int LOGN>
double sup_norm(const Message<LOGN> &z){
    double max = 0;
    for(int i =0 ; i < (1 << (LOGN - 1)); ++i) {
        double abs = sqrt(z.r[i] * z.r[i] + z.i[i] * z.i[i]);
        max = abs > max ? abs : max;
    }
    return max;
}

template <int LOGN>
void print(const std::string name, const Message<LOGN> &z){
	std::cout << "Message " << name << std::endl;
	for(int i=0; i <min((1 << (LOGN - 1)), 10); ++i) {
		std::cout << "(" << z.r[i] << ", " << z.i[i] << ") ";
	}
	std::cout << std::endl;
}

template <int LOGN>
void conj(const Message<LOGN> &z1, Message<LOGN> &z2) {
    for(int i = 0; i < (1 << (LOGN - 1)); ++i) {
        z2.r[i] = z1.r[i];
        z2.i[i] = -z1.i[i];
    }
}

template<int LOGN, int LOGDELTA, int G = 1>
void CoeffToSlot (const Message<LOGN>& z,
					Message<LOGN> z_[2]){
	const int N = 1 << LOGN;
	static bool is_init = false;
	static SparseDiagonal<N/2,3> U0r[LOGN-1];
    static SparseDiagonal<N/2,3> U0i[LOGN-1];

	if(!is_init) {
		splitU0NR<LOGN>(U0r, U0i);
    	for (int n = 0; n < LOGN - 1; n++) {
			U0r[n].transpose();
			U0i[n].transpose();
			U0i[n].negate();
			U0r[n] *= 1 / pow(N, 1.0 / (LOGN - 1));
			U0i[n] *= 1 / pow(N, 1.0 / (LOGN - 1));
		}
		is_init = true;
	}

    Message <LOGN> z1, U0z;
    z1 = z;
    for(int d = 0; d < LOGN - 1; ++d) {
		U0z = z1;
        matrix_vector_product<LOGN, 3>(U0z, U0r[d], U0i[d], z1);
	}

    Message<LOGN> zi;
    for(int i = 0; i < (1 << (LOGN - 1)); ++i) {
        zi.r[i] = 0;
        zi.i[i] = 1;
    }

    Message <LOGN> z2;
    conj(z1, z2);
    add(z1, z2, z_[0]);
    sub(z1, z2, z_[1]);
    mult(z_[1], zi, z_[1]);
}