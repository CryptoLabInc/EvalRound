#pragma once

#include "util.h"

#include <cstdint>
#include <cmath>

constexpr int LOGDELTA = 10;
constexpr std::uint64_t Delta = 1ULL << LOGDELTA;
constexpr std::uint64_t DeltaSq = 1ULL << (2 * LOGDELTA);
constexpr std::uint64_t DeltaTr = 1ULL << (3 * LOGDELTA);
constexpr int LOGQ = 24;
constexpr int LOGQ_UP = LOGQ * 2;
constexpr int N = 1 << 5;
constexpr int h = 64;
constexpr int K = 16;

void set_test_message(double zr[N/2], double zi[N/2]) {
    sampleUniform(zr, zr+N/2);
    sampleUniform(zi, zi+N/2);
/*    for(int i = 0; i < N / 2; ++i) {
        double x = (double) (i*i) / (N/2) * M_PI;
        zr[i] = cos(x);
        zi[i] = sin(x);
    }
    */
}

void set_test_matrix(double Ar[K][N/2], double Ai[K][N/2]) {
    /*for(int k = 0; k < K; ++k) {
        for(int i = 0; i < N/2; ++i) {
            double x = (double) (i + k) / (N/2) * M_PI;
            Ar[k][i] = cos(x);
        }
    }

    for(int k = 0; k < K; ++k) {
        for(int i = 0; i < N/2; ++i) {
            double x = (double) (i + k) / (N/2) * M_PI;
            Ai[k][i] = sin(x);
        }
    }*/

    for(int k = 0; k < K; ++k) {
        sampleUniform(Ar[k], Ar[k]+N/2);
        sampleUniform(Ai[k], Ai[k]+N/2);
    }
}