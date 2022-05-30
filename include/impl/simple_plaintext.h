#pragma once

#include "message.h"
#include "HEAAN/matrix/DFT.h"
#include "util/util.h"

#include <iostream>
#include <string>

template <int LOGN>
struct SimplePlaintext {
    double data[1 << LOGN] = {0};

    double operator[]( int i ) const{ return data[i]; }
    double& operator[]( int i ) { return data[i]; } 
};

template <int LOGN>
void encode(const Message<LOGN> &z, uint64_t Delta, SimplePlaintext<LOGN> &pt){
	const int N = 1 << LOGN;
	double m[N]; ifft<LOGN>(z.r,z.i,m);
	double Delta_double = (double) Delta;
	for(int i=0; i<N; i++){
		pt[i] = (double) round(m[i]*Delta_double);
	}
}

template <int LOGN>
void encode_raw(const Message<LOGN> &z, uint64_t Delta, SimplePlaintext<LOGN> &pt){
	const int N = 1 << LOGN;
	double m[N]; ifft<LOGN>(z.r,z.i,m);
	double Delta_double = (double) Delta;
	for(int i=0; i<N; i++){
		pt[i] = m[i]*Delta_double;
	}
}

template <int LOGN>
void decode(const SimplePlaintext<LOGN> &pt, uint64_t Delta, Message<LOGN> &z){
	const int N = 1 << LOGN;
	double m[N];
	double Delta_double = (double) Delta;
    for(int i = 0; i < N; ++i) {
        m[i] = pt[i] / Delta_double;
    }
	fft<LOGN>(m,z.r,z.i);
}

template <int LOGN>
void conv(const SimplePlaintext<LOGN> &pt1, const SimplePlaintext<LOGN> &pt2, SimplePlaintext<LOGN> &pt3){
	    const int N = 1 << LOGN;
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int k=0;   k<=i; k++) sum += pt1[k]*pt2[i-k];
		for(int k=i+1; k< N; k++) sum -= pt1[k]*pt2[i+N-k];
		pt3[i]=sum;
	}
}

template <int LOGN>
void add(const SimplePlaintext<LOGN> &pt1, const SimplePlaintext<LOGN> &pt2, SimplePlaintext<LOGN> &pt3){
	for(int i=0; i<(1 << LOGN); i++){
		pt3[i] = pt1[i] + pt2[i];
	}
}

template <int LOGN>
void sub(const SimplePlaintext<LOGN> &pt1, const SimplePlaintext<LOGN> &pt2, SimplePlaintext<LOGN> &pt3){
	for(int i=0; i<(1 << LOGN); i++){
		pt3[i] = pt1[i] - pt2[i];
	}
}

template <int LOGN>
double square_sum(const SimplePlaintext<LOGN> &pt) {
	double sum = 0;
	for(int i = 0; i < (1 << LOGN); ++i) {
		sum += pt[i]* pt[i];
	}
	return sum;
}

template <int LOGN>
double norm(const SimplePlaintext<LOGN> &pt){
    return sqrt(square_sum(pt));
}

template <int LOGN>
void rotate_once(const SimplePlaintext<LOGN> &pt, SimplePlaintext<LOGN> &pt_rot){
	const int N = 1 << LOGN;
	for(int i=0; i<N; i++){
		int q = (5*i)/N;
		int r = (5*i)%N;
		pt_rot[r] = pt[i];
		if(q%2==1) pt_rot[r] = -pt_rot[r];
	}
}

template <int LOGN>
void rotate(const SimplePlaintext<LOGN> &pt, SimplePlaintext<LOGN> &pt_rot, int r){
	const int N = 1 << LOGN;
	SimplePlaintext<LOGN> pt_tmp;
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

template <int LOGN>
void print(const std::string name, const SimplePlaintext<LOGN> &pt){
	std::cout << "SimplePlaintext " << name << std::endl;
	for(int i=0; i <min(1 << LOGN, 10); ++i) {
		std::cout << pt[i] << " ";
	}
	std::cout << std::endl;
}