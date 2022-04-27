#pragma once

#include "define.h"

#include <iostream>
#include <string>


void rotate_message(const double z[N/2], double z_rot[N/2], int r) {
    for(int i = 0; i < N/2; ++i)
        z_rot[i] = z[(i + r) % (N/2)];
}

void matrix_vector_product(
    const double zr[N/2], const double zi[N/2],
    const double Ar[K][N/2], const double Ai[K][N/2],
    double Azr[N/2], double Azi[N/2]) {
    for(int i = 0; i < N/2; ++i) {
        double sumr = 0, sumi = 0;
        for(int k = 0; k < K; ++k) {
            int j = (i+k) % (N/2);
            sumr += Ar[k][i] * zr[j] - Ai[k][i] * zi[j];
            sumi += Ar[k][i] * zi[j] + Ai[k][i] * zr[j];
        }
        Azr[i] = sumr;
        Azi[i] = sumi;
    }
}

double square_sum(const double zr[N/2], const double zi[N/2]){
    double sum = 0;
    for(int i = 0; i < N/2; ++i) {
        sum += zr[i] * zr[i] + zi[i] * zi[i];
    }
    return sum;
}

double norm(const double zr[N/2], const double zi[N/2]){
    return sqrt(square_sum(zr, zi));
}

void print(const std::string name , const double z[N/2]){
	std::cout << "Message " << name << std::endl;
	for(int i=0; i <N/2; ++i) {
		std::cout << z[i] << " ";
	}
	std::cout << std::endl;
}