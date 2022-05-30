#pragma once

#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <set>

#include <chrono>
#include <string>

#define PI 3.1415926535897932384626433

inline int log(int N) {
    int i = 0, n = 1;
    while(n < N) {
        i++;
        n *= 2;
    }
    return i;
}

inline uint32_t bitReverse32(uint32_t x) {
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return ((x >> 16) | (x << 16));
}

inline uint32_t bitReverse(uint32_t x, int digits) {
    return bitReverse32(x) >> (32 - digits);
}

uint64_t rand_uint64(void){
	uint64_t r =0;
	for( int i=0; i<64; i+=15 /*30*/){
		  r = r*((uint64_t)RAND_MAX + 1) + rand();
  }
  return r;
}

template<class T>
T min(T a, T b){
    return a < b ? a : b;
}

template<int N>
void print_array(const uint64_t A[N]) {
    for(int i = 0; i < min(N, 10); ++i)
        std::cout << A[i] << " ";
    std::cout << std::endl;
}

template<int N>
void print_array(const double A[N]) {
    for(int i = 0; i < min(N, 10); ++i)
        std::cout << A[i] << " ";
    std::cout << std::endl;
}

#include <chrono>

class Timer {
    public:
        Timer(std::string name) {
            name_ = name;
            begin_ = std::chrono::steady_clock::now();
            std::cout << name_ << " start" << std::endl;
        }

        void stop() {
            end_ = std::chrono::steady_clock::now();
            std::cout << name_ << " : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_ - begin_).count() << "[ms]" << std::endl;
        }

    private:
        std::string name_;
        std::chrono::steady_clock::time_point begin_;
        std::chrono::steady_clock::time_point end_;
};
