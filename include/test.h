#pragma once

#include "define.h"
#include "message.h"
#include "plaintext.h"

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

void set_test_rounded_message(double zr[N/2], double zi[N/2]) {
    set_test_message(zr, zi);
    double pt[N];
    encode(zr, zi, Delta, pt);
    decode(pt, Delta, zr, zi);
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