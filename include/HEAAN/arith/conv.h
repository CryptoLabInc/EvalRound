#pragma once

#include "conv/conv.h"

template<int N>
void conv( const int s1[N],
		   const int s2[N], int s3[N]){
	for(int i=0; i<N; i++){
		s3[i]=0;
		for(int k=0;   k<=i; k++) s3[i] += s1[k]*s2[i-k];
		for(int k=i+1; k< N; k++) s3[i] -= s1[k]*s2[i+N-k];
	}
}

template< int LOGQ, int N >
void conv( const int s[N],
	       const R_Q<LOGQ,N>& A, R_Q<LOGQ,N>& C ){
	#pragma omp parallel for
	for(int i=0;i<N;i++){
		C[i].setzero(); const Z_Q<LOGQ>* temp;
		for(int k=0;k<N;k++){
			if( s[k]!=0 ){
				bool sign=s[k]==1;
				if(k<=i) temp=&A[  i-k]; 
				else    {temp=&A[N+i-k]; sign=!sign;}
			
				if(sign) C[i]+=*temp;
				else     C[i]-=*temp;
			}
		}
	}
}

template< int LOGQ, int N >
void conv( const R_Q<LOGQ,N>& A,
	       const R_Q<LOGQ,N>& B, R_Q<LOGQ,N>& C ){
    // we fix params for primes since this will be used only for conv.
    // 102 60bit primes, mod(q[i], 2^18) = 1
    static constexpr uint64_t Q_primes[] = {
        0xffffffffffc0001ULL,
        0xfffffffff840001ULL,
        0x1000000000980001ULL,
        0x1000000000b00001ULL,
        0xfffffffff240001ULL,
        0x1000000000f00001ULL,
        0xffffffffe7c0001ULL,
        0xffffffffe740001ULL,
        0x1000000001a00001ULL,
        0xffffffffe4c0001ULL,
        0xffffffffe440001ULL,
        0xffffffffe400001ULL,
        0x1000000002340001ULL,
        0xffffffffdbc0001ULL,
        0xffffffffd840001ULL,
        0x1000000002940001ULL,
        0xffffffffd680001ULL,
        0xffffffffd000001ULL,
        0xffffffffcf00001ULL,
        0xffffffffcdc0001ULL,
        0xffffffffcc40001ULL,
        0x1000000003680001ULL,
        0x1000000003900001ULL,
        0xffffffffc300001ULL,
        0x1000000003ec0001ULL,
        0xffffffffbf40001ULL,
        0xffffffffbdc0001ULL,
        0xffffffffb880001ULL,
        0x1000000004d40001ULL,
        0x1000000004f80001ULL,
        0x1000000005040001ULL,
        0xffffffffaec0001ULL,
        0x1000000005400001ULL,
        0x1000000005480001ULL,
        0x10000000054c0001ULL,
        0x1000000005640001ULL,
        0xffffffffa380001ULL,
        0xffffffffa200001ULL,
        0xffffffffa0c0001ULL,
        0x10000000069c0001ULL,
        0xffffffff9600001ULL,
        0x1000000006e00001ULL,
        0xffffffff91c0001ULL,
        0xffffffff8f40001ULL,
        0x10000000075c0001ULL,
        0x1000000007640001ULL,
        0xffffffff8680001ULL,
        0x1000000007d40001ULL,
        0x10000000080c0001ULL,
        0xffffffff7e40001ULL,
        0xffffffff7bc0001ULL,
        0x1000000008540001ULL,
        0xffffffff76c0001ULL,
        0xffffffff7680001ULL,
        0x1000000008e40001ULL,
        0xffffffff6fc0001ULL,
        0x1000000009240001ULL,
        0x1000000009300001ULL,
        0xffffffff6880001ULL,
        0x10000000098c0001ULL,
        0x1000000009c00001ULL,
        0xffffffff6340001ULL,
        0xffffffff5d40001ULL,
        0x100000000a8c0001ULL,
        0xffffffff54c0001ULL,
        0x100000000ab80001ULL,
        0x100000000aec0001ULL,
        0x100000000af40001ULL,
        0xffffffff4d40001ULL,
        0x100000000b300001ULL,
        0x100000000b400001ULL,
        0x100000000b7c0001ULL,
        0x100000000ba00001ULL,
        0xffffffff4380001ULL,
        0xffffffff3e80001ULL,
        0xffffffff37c0001ULL,
        0xffffffff36c0001ULL,
        0x100000000d400001ULL,
        0x100000000d4c0001ULL,
        0x100000000db00001ULL,
        0x100000000de00001ULL,
        0xffffffff2100001ULL,
        0x100000000df40001ULL,
        0x100000000e1c0001ULL,
        0xffffffff1d80001ULL,
        0xffffffff1cc0001ULL,
        0x100000000e3c0001ULL,
        0x100000000e480001ULL,
        0xffffffff1900001ULL,
        0xffffffff1740001ULL,
        0xffffffff15c0001ULL,
        0x100000000f140001ULL,
        0xffffffff0e80001ULL,
        0x100000000f200001ULL,
        0x100000000f3c0001ULL,
        0xfffffffeff80001ULL,
        0xfffffffeff40001ULL,
        0x1000000010100001ULL,
        0x1000000010a00001ULL,
        0xfffffffeefc0001ULL,
        0x1000000011340001ULL,
        0xfffffffee8c0001ULL
    };
    const int L = (2*LOGQ + 17 + 59) / 60;
	static CONV<L, LOGQ, N> *theConv = nullptr;
    static bool theConv_initialized = false;
    if(theConv_initialized == false) {
        theConv = new CONV<L, LOGQ, N>(Q_primes);
        theConv_initialized = true;
    }
    theConv->conv(A, B, C);
}