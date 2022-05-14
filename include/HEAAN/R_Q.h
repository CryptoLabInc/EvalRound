#pragma once                
#include "Z_Q.h"            

//----------------------------------------------
// This is an implementation of the ring (+,-,*)
// Z_Q[x]/x^N+1
//----------------------------------------------
template <int LOGQ, int N>
struct R_Q{
	   Z_Q<LOGQ> coeff[N];

	   R_Q() {}
	   R_Q           (const R_Q<LOGQ, N>& A){for(int i=0;i<N;i++) coeff[i] = A.coeff[i]; }
	   void operator=(const R_Q<LOGQ, N>& A){for(int i=0;i<N;i++) coeff[i] = A.coeff[i]; }
	   const Z_Q<LOGQ>& operator[](int i) const {return coeff[i];}
	         Z_Q<LOGQ>& operator[](int i)       {return coeff[i];}
	   void operator += (const R_Q<LOGQ,N>& B);
	   void operator -= (const R_Q<LOGQ,N>& B);
	   void operator *= ( uint64_t a );
	   void setzero();
	   void print_unsigned( )const;
	   void print_signed  ( )const;
};
//-----------------------------------------------------------
// negative-wrapped convolution
//-----------------------------------------------------------

template< int LOGQ, int N >
void conv( const int s[N],
	       const R_Q<LOGQ,N>& A, R_Q<LOGQ,N>& C ){
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
//-----------------------------------------------------------
//
//
//-----------------------------------------------------------
template <int LOGQ, int N>
void R_Q <LOGQ,N> :: operator += (const R_Q<LOGQ,N>& B){
	for (int i=0; i<N; i++)
		coeff[i]+=B.coeff[i];
}
template <int LOGQ, int N>
void R_Q<LOGQ, N> :: operator -= (const R_Q<LOGQ, N>& B){
	for (int i=0; i<N; i++)
		coeff[i]-=B.coeff[i];
}

template <int LOGQ, int N>
void R_Q <LOGQ,N> :: operator *= (uint64_t a){
	for (int i=0; i<N; i++)
		coeff[i]*=a;
}

template <int LOGQ, int N>
void R_Q<LOGQ, N>::print_unsigned()const{
	printf("[");
	for(int i=0;i<N;i++){
		coeff[i].print_unsigned();
		if(i<N-1) printf(",");
	}
	printf("]\n");
}

template <int LOGQ, int N>
void R_Q<LOGQ, N>::print_signed()const{
	printf("[");
	for(int i=0;i<N;i++){
		if(coeff[i].is_bigger_than_halfQ()){
			printf("-"); Z_Q<LOGQ> temp=coeff[i];
			temp.negate();
			temp.print_unsigned();
		}
		else
			coeff[i].print_unsigned();
		if(i<N-1) printf(",");
	}
	printf("]\n");
}

template <int LOGQ, int N>
void R_Q<LOGQ, N>::setzero(){
	for(int i=0;i<N;i++)
		coeff[i].setzero();
}