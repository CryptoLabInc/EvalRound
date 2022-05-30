#pragma once                
#include "R_Q.h"  
#include "conv.h"          
//----------------------------------------------
// This is an implementation of the ring (+,-,*)
// R_Q^2[x]/x^N+1
//----------------------------------------------
template <int LOGQ, int N>
struct R_Q_square{
	   R_Q<LOGQ,N> data[2];
	   R_Q_square() {}
	   R_Q_square    (const R_Q_square<LOGQ,N>& A){data[0]=A.data[0]; data[1]=A.data[1]; }
	   void operator=(const R_Q_square<LOGQ,N>& A){data[0]=A.data[0]; data[1]=A.data[1]; }
	   const R_Q<LOGQ,N>& operator[](int i)const {return data[i];}
			 R_Q<LOGQ,N>& operator[](int i)		 {return data[i];}
	   void operator += (const R_Q_square<LOGQ,N>& B);
	   void operator -= (const R_Q_square<LOGQ,N>& B);
	   void operator *= (const R_Q       <LOGQ,N>& B);
	   void operator *= (uint64_t a );
	   void print (    )const;
	   void setzero(){ data[0].setzero(); data[1].setzero();}
};
//-----------------------------------------------------------
//  ���
//-----------------------------------------------------------
template <int LOGQ, int N>
void R_Q_square <LOGQ, N>::print()const{
	printf("[");
	data[0].print_unsigned(); printf(",\n");
	data[1].print_unsigned();
	printf("]\n");
}
//-----------------------------------------------------------
// �⺻ operations
//-----------------------------------------------------------
template<int LOGQ, int N>
void R_Q_square <LOGQ,N> :: operator += (const R_Q_square<LOGQ,N>& B){
	data[0]+=B.data[0];
	data[1]+=B.data[1];
}

template<int LOGQ, int N>
void R_Q_square <LOGQ,N> :: operator -= (const R_Q_square<LOGQ,N>& B){
	data[0]-=B.data[0];
	data[1]-=B.data[1];
}


template<int LOGQ, int N>
void R_Q_square <LOGQ,N> :: operator *= (const R_Q<LOGQ,N>& B){
	R_Q<LOGQ,N> temp;
	temp=data[0]; conv<LOGQ,N>(temp,B,data[0]);
	temp=data[1]; conv<LOGQ,N>(temp,B,data[1]);
}


template<int LOGQ, int N>
void R_Q_square <LOGQ,N> :: operator *= (uint64_t a){
	data[0]*=a;
	data[1]*=a;
}

template<int LOGQfrom, int LOGQto, int N>
void resize(const R_Q_square<LOGQfrom, N> &Afrom, R_Q_square<LOGQto, N> &Ato) {
	for(int i = 0; i < 2; ++i) {
		resize(Afrom.data[i], Ato.data[i]);
	}
}

