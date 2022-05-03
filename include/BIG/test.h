#pragma once

#include "define.h"
#include "message.h"
#include "plaintext.h"
#include "util.h"

#include "HEAAN/DFT.h"

void set_test_message(double zr[N/2], double zi[N/2]) {
    //sampleUniform(zr, zr+N/2);
    //sampleUniform(zi, zi+N/2);
    for(int i = 0; i < N / 2; ++i) {
        double x = (double) (i) / (N/2) * M_PI;
        zr[i] = cos(x) / sqrt(2);
        zi[i] = sin(x) / sqrt(2);
    }
}

void set_random_message(double zr[N/2], double zi[N/2]) {
    sampleUniform(zr, zr+N/2);
    sampleUniform(zi, zi+N/2);
}

void set_test_message(double m[N]) {
    //sampleUniform(m, m+N);
    for(int i = 0; i < N; ++i) {
        double x = (double) (i) / N * M_PI;
        m[i] = cos(x);
    }
    
}

void set_test_rounded_message(double zr[N/2], double zi[N/2]) {
    set_test_message(zr, zi);
    R_Q<LOGQ, N> pt;
    encode<LOGQ, LOGN>(zr, zi, Delta, pt);
    decode<LOGQ, LOGN>(pt, Delta, zr, zi);
}

void set_test_matrix(double Ar[K][N/2], double Ai[K][N/2]) {
    for(int k = 0; k < K; ++k) {
        //sampleUniform(Ar[k], Ar[k]+N/2);
        //sampleUniform(Ai[k], Ai[k]+N/2);
        for(int i = 0; i < N / 2; ++i) {
            double x = (double) (k*i) / (N/2) * M_PI;
            Ar[k][i] = cos(x) / sqrt(2);
            Ai[k][i] = sin(x) / sqrt(2);
        }
    }
}

void set_random_matrix(double Ar[K][N/2], double Ai[K][N/2]) {
    for(int k = 0; k < K; ++k) {
        sampleUniform(Ar[k], Ar[k]+N/2);
        sampleUniform(Ai[k], Ai[k]+N/2);
    }
}

void set_test_U0_matrix(SparseDiagonal<(1<<(LOGN-1)),3> U0r[LOGN-1],
	            SparseDiagonal<(1<<(LOGN-1)),3> U0i[LOGN-1]) {
    splitU0NR<LOGN>(U0r, U0i);
    
    for (int n = 0; n < LOGN -1; n++) {
		U0r[n].transpose();
		U0i[n].transpose();
		U0i[n].negate();
	}
}
