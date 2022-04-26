#pragma once

#include "HEAAN/matrix.h"

#include <cstdint>
#include <cmath>

constexpr int LOGDELTA = 3;
constexpr std::uint64_t Delta = 1ULL << LOGDELTA;
constexpr std::uint64_t DeltaSq = 1ULL << (2 * LOGDELTA);
constexpr int LOGQ = 24;
constexpr int LOGQ_UP = LOGQ * 2;
constexpr int N = 1 << 5;
constexpr int h = 64;
constexpr int K = 1;

void set_test_message(double zr[N/2], double zi[N/2]) {
    for(int i = 0; i < N / 2; ++i) {
        double x = (double) i / (N/2) * M_PI;
        zr[i] = cos(x);
        zi[i] = sin(x);
    }
}

void set_test_matrix(SparseDiagonal<N/2, K> &Ar, SparseDiagonal<N/2, K> &Ai) {
    for(int k = 0; k < K; ++k) {
        Ar.off[k] = k;
        Ar.zero[k] = false;
        for(int i = 0; i < N/2; ++i) {
            double x = (double) (i + k) / (N/2) * M_PI;
            Ar.vec[k][i] = cos(x);
        }
    }

    for(int k = 0; k < K; ++k) {
        Ai.off[k] = k;
        Ai.zero[k] = false;
        for(int i = 0; i < N/2; ++i) {
            double x = (double) (i + k) / (N/2) * M_PI;
            Ai.vec[k][i] = sin(x);
        }
    }
}

void rotate_message(const double z[N/2], double z_rot[N/2], int r) {
    for(int i = 0; i < N/2; ++i)
        z_rot[i] = z[(i + r) % (N/2)];
}

void matrix_vector_product(
    const double zr[N/2], const double zi[N/2],
    const SparseDiagonal<N/2, K> Ar, const SparseDiagonal<N/2, K> Ai,
    double Azr[N/2], double Azi[N/2]) {
    for(int i = 0; i < N/2; ++i) {
        double sumr = 0, sumi = 0;
        for(int k = 0; k < K; ++k) {
            sumr += Ar.vec[k][i] * zr[i] - Ai.vec[k][i] * zi[i];
            sumi += Ar.vec[k][i] * zi[i] + Ai.vec[k][i] * zr[i];
        }
        Azr[i] = sumr;
        Azi[i] = sumi;
    }
}