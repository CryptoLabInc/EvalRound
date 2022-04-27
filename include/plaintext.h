#pragma once

#include <iostream>
#include <string>

#include "HEAAN/DFT.h"
#include "define.h"

void encode( const double zr[N/2], 
			 const double zi[N/2], uint64_t Delta, int64_t pt[N]){
	double m[N]; idft<N>(zr,zi,m);
	double Delta_double = (double) Delta;
	for(int i=0; i<N; i++){
		pt[i] = (int64_t) round(m[i]*Delta_double);
	}
}

void decode( const int64_t pt[N], uint64_t Delta, double zr[N/2], 
													double zi[N/2]){
	double m[N];
	double Delta_double = (double) Delta;
    for(int i = 0; i < N; ++i) {
        m[i] = pt[i] / Delta_double;
    }
	dft<N>(m,zr,zi);
}

void conv(const int64_t pt1[N], const int64_t pt2[N], int64_t pt3[N]){
	for(int i=0; i<N; i++){
		int64_t sum = 0;
		for(int k=0;   k<=i; k++) sum += pt1[k]*pt2[i-k];
		for(int k=i+1; k< N; k++) sum -= pt1[k]*pt2[i+N-k];
		pt3[i]=sum;
	}
}

void set_zero(int64_t pt[N]) {
	for(int i = 0; i <N; ++i) {
		pt[i] = 0;
	}
}

void add(const int64_t pt1[N], const int64_t pt2[N], int64_t pt3[N]){
	for(int i=0; i<N; i++){
		pt3[i] = pt1[i] + pt2[i];
	}
}

void sub(const int64_t pt1[N], const int64_t pt2[N], int64_t pt3[N]){
	for(int i=0; i<N; i++){
		pt3[i] = pt1[i] - pt2[i];
	}
}

double square_sum(const int64_t pt[N]) {
	double sum = 0;
	for(int i = 0; i < N; ++i) {
		sum += (double) pt[i]* (double) pt[i];
	}
	return sum;
}

double norm(const int64_t pt[N]){
    return sqrt(square_sum(pt));
}

void print(const std::string name, const int64_t pt[N]){
	std::cout << "Plaintext " << name << std::endl;
	for(int i=0; i <N; ++i) {
		std::cout << pt[i] << " ";
	}
	std::cout << std::endl;
}

void rotate_once(const int64_t pt[N], int64_t pt_rot[N]){
	for(int i=0; i<N; i++){
		int q = (5*i)/N;
		int r = (5*i)%N;
		pt_rot[r] = pt[i];
		if(q%2==1) pt_rot[r] = -pt_rot[r];
	}
}

void rotate(const int64_t pt[N], int64_t pt_rot[N], int r){
	int64_t pt_tmp[N];
	for(int i = 0; i < N; ++i) {
		pt_rot[i] = pt[i];
	}
	for(int rot = 0; rot < r; ++rot) {
		for(int i = 0; i < N; ++i) {
			pt_tmp[i] = pt_rot[i];
		}
		rotate_once(pt_tmp, pt_rot);
	}
}

void matrix_vector_product(
    const int64_t pt[N],
    const double Ar[K][N/2], const double Ai[K][N/2],
    int64_t pt_Az[N/2]){
	int64_t pt_v[N], pt_rot[N], pt_conv[N];
	set_zero(pt_Az);
	for(int k = 0; k < K; ++k) {
		encode(Ar[k], Ai[k], Delta, pt_v);
		rotate(pt, pt_rot, k);
		conv(pt_v, pt_rot, pt_conv);
		add(pt_Az, pt_conv, pt_Az);
	}
}