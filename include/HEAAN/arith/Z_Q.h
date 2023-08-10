#pragma once

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "util/util.h"

//--------------------------------------------------------------------
//
// A simple implemenation of the ring Z_Q
//
// Usage : Z_Q<1270> A, B;
//         A+=B;
//         A.print_unsigned();
//
//--------------------------------------------------------------------
template< int LOGQ >
struct Z_Q {
	//-----------------------------------------------------------
	// Q=2^LOGQ, memory functions
	//-----------------------------------------------------------
	uint64_t data[ (LOGQ+63)/64];
	Z_Q(){}
	Z_Q           ( const Z_Q<LOGQ>& ); // copy constructor
	void operator=( const Z_Q<LOGQ>& ); // copy operator

	uint64_t  operator[]( int i ) const{ return data[i]; }
	uint64_t& operator[]( int i )      { return data[i]; }
	int get_length() const{ return (LOGQ+63)/64; }
	void fr_str( const char* str )     ; // string of decimal digits
	void to_str(       char* str )const; //    str of size>=520

	//-----------------------------------------------------------
	// basic operations
	//-----------------------------------------------------------
	void operator+=(const Z_Q<LOGQ>& A); 
	void operator-=(const Z_Q<LOGQ>& A);
	void operator*=(const Z_Q<LOGQ>& A);
	void operator*=(uint64_t a);
	void remove_clutter();
	void negate();
	void setzero();

	//-----------------------------------------------------------
	// etc
	//-----------------------------------------------------------
	void print_unsigned()const;
	bool is_bigger_than_halfQ()const;
	int max_valid_digit() const;
	explicit operator double() const;
};

//-----------------------------------------------------------
// shift
//-----------------------------------------------------------
template< int LOGQfr, int LOGQto>
void shift_right( const Z_Q<LOGQfr>& A,
	                    Z_Q<LOGQto>& B ){
	if(LOGQfr == LOGQto) {
		for (int i = 0; i < (LOGQto + 63) / 64; ++i)
			B[i] = A[i];
		return;
	}
	int q=(LOGQfr-LOGQto)/64;
	int r=(LOGQfr-LOGQto)%64;
	for(int i=0;i<(LOGQto+63)/64;i++)
		B[i]=A[i+q];
	if(r > 0) {
		uint64_t Bextra=0;
		if((LOGQto+r+63)/64 != (LOGQto+63)/64)
			Bextra = A[(LOGQto+63)/64+q];
		Bextra<<=64-r;
		for(int i=(LOGQto+63)/64-1;i>=0;i--){
			uint64_t temp=B[i]<<(64-r);
			B[i]>>=r;
			B[i]+=Bextra; Bextra=temp;
		}
	}
	int idx = LOGQfr - LOGQto - 1;
	bool c = ((uint64_t)0x1) & (A[idx / 64] >> (idx % 64));
	if(c) {
		Z_Q<LOGQto> C;
		for (int i = 0; i < (LOGQto + 63) / 64; ++i)
			C[i] = 0;
		C[0] = 1;
		B += C;
		return;
	}
}

template< int LOGQfr, int LOGQto>
void shift_left ( const Z_Q<LOGQfr>& A,
	                    Z_Q<LOGQto>& B ){
	int q=(LOGQto-LOGQfr)/64;
	int r=(LOGQto-LOGQfr)%64;
	for(int i=0;i<(LOGQfr+63)/64;i++)
		B[i+q]=A[i];
	for(int i=0;i<q;i++) B[i]=0;
	if(r == 0)
		return;
	uint64_t Bextra=0;
	for(int i=0;i<(LOGQfr+63)/64;i++){
		uint64_t temp=B[i+q]>>(64-r);
		B[i+q]<<=r;
		B[i+q]+=Bextra; Bextra=temp;
	}
	if((LOGQfr+r+63)/64 != (LOGQfr+63)/64)
		B[(LOGQfr+63)/64+q]=Bextra;
}
//-----------------------------------------------------------
// implementations of memory functions
//-----------------------------------------------------------
template< int LOGQ>
Z_Q<LOGQ>::Z_Q( const Z_Q<LOGQ>& B){ // copy constructor
	for(int i=0;i<get_length();i++) data[i]=B.data[i];
}

template< int LOGQ>
void Z_Q<LOGQ>::operator=( const Z_Q<LOGQ>& B){ // copy operator
	for(int i=0;i<get_length();i++) data[i]=B.data[i];
}

template< int LOGQ >
void Z_Q<LOGQ>::fr_str( const char* str ){ // string of decimal digits
	// read decimal digits in the reverse direction
	int A[520];
	int n=strlen(str);
	for(int i=0;i<n;i++)
		A[i]=str[n-1-i]-'0';
	// from base 10 to base 2^32
	uint64_t B[((LOGQ+63)/64)*2];
	bool zero=false; int count=0;
	while(!zero){
		uint64_t rem=0; zero=true;
		for(int i=n-1;i>=0;i--){
			rem*=10; rem+=A[i];
			A[i]=rem>>32; if(A[i]!=0) zero=false;
			rem = (rem<<32)>>32;
		}
		B[count]=rem; count++;
	}
	for(int i=count;i<2*get_length();i++) B[i]=0;
	// from base 2^32 to base 2^64
	for(int i=0;i<get_length();i++)
		data[i] = B[2*i]+(B[2*i+1]<<32);
}

template< int LOGQ >                        // str of size>=520
void Z_Q<LOGQ>::to_str( char* str ) const { // string of decimal digits
	const int length = (LOGQ+63)/64;
	// base 2^64 => 2^32 
	uint64_t A[2*length];
	for(int i=0;i<length;i++){
		A[2*i  ] = (data[i]<<32)>>32;
		A[2*i+1] = (data[i]>>32)    ;
	}
	// base 2^32 => base 10
	assert(LOGQ <= 1300 );
	int B[520]; int count=0;
	bool zero=false;
	while(!zero){
		uint64_t rem=0; zero=true;
		for(int i=2*length-1;i>=0;i--){
			rem<<=32; rem+=A[i];
			A[i]=rem/10; if(A[i]!=0) zero=false;
			rem =rem%10;
		}
		B[count]=(int)rem; count++;
	}
	for(int i=0;i<count;i++) str[i]='0'+B[count-1-i];
	str[count]=0;
}


//-----------------------------------------------------------
// ������ ���
//-----------------------------------------------------------
template<int LOGQ>
void Z_Q<LOGQ>::print_unsigned()const {
	char str[520]; to_str(str);
	printf("%s\n",str);
}

//-----------------------------------------------------------
// �⺻ operations
//-----------------------------------------------------------
template<int LOGQ>
void Z_Q<LOGQ>::operator+=(const Z_Q<LOGQ>& A) {
	int length = (LOGQ + 63) / 64; bool c = 0;
	for (int i = 0; i < length; i++) {
		data[i] += A.data[i];
		bool c_next = data[i] < A.data[i];
		if (c) data[i]++;
		c = (c&&data[i] == 0) || c_next;
	}
	remove_clutter();
}

template<int LOGQ>
void  Z_Q<LOGQ>::remove_clutter() {
	int length = (LOGQ + 63) / 64;
	int r = 64 * length - LOGQ;
	data[length - 1] <<= r;
	data[length - 1] >>= r;
}


template<int LOGQ>
void  Z_Q<LOGQ>::negate() {
	int length = (LOGQ + 63) / 64;
	uint64_t ones = -1;
	for (int i = 0; i < length; i++) data[i] = data[i] ^ ones;
	remove_clutter();
	data[0]++; bool c = (data[0] == 0);
	for (int i = 1; c && (i < length); i++) {
		data[i]++;
		c = (data[i] == 0);
	}
	remove_clutter();
}

template<int LOGQ>
void Z_Q<LOGQ>::setzero(){
	int length = (LOGQ + 63) / 64;
	for(int i=0;i<length;i++) data[i]=0;
}

template<int LOGQ>
void  Z_Q<LOGQ>::operator-=(const Z_Q<LOGQ>& A) {
	Z_Q<LOGQ> negA(A);
	negA.negate();
	*this+=(negA);
}

void mul(uint64_t a, uint64_t b, uint64_t& lo, uint64_t& hi){
	uint64_t ah = a>>32, al = (a<<32)>>32;
	uint64_t bh = b>>32, bl = (b<<32)>>32;
	uint64_t cll =al*bl;
	uint64_t clh =al*bh;
	uint64_t chl =ah*bl;
	uint64_t chh =ah*bh;

	chl += (cll>>32)+((clh<<32)>>32);
	chh += (chl>>32)+(clh>>32);
	hi=chh;
	lo=((cll<<32)>>32)+(chl<<32);
}

template<int LOGQ>
void Z_Q<LOGQ>::operator*=(uint64_t a){
	Z_Q<LOGQ>Hi; Hi.data[0]=0;
	for(int i=0; i<get_length(); i++){
		if(i<get_length()-1)
			mul(data[i],a,   data[i  ],
						  Hi.data[i+1]);
		else{ uint64_t temp;
			mul(data[i],a,   data[i  ],
						  temp);
		}
	}
	*this += Hi;
}

template<int LOGQ>
void  Z_Q<LOGQ>::operator*=(const Z_Q<LOGQ>& B) {
	Z_Q<LOGQ> sum; sum.setzero();
	Z_Q<LOGQ> ABi; int length = get_length();
	for(int i=0; i<length; i++){
		ABi = *this;
		ABi *= B.data[i];
		for(int j=length-1; j>=0; j--)
			if(j>=i) ABi.data[j]=ABi.data[j-i];
			else	 ABi.data[j]=0;
		
		sum += ABi;
	}
	*this=sum;
}


template<int LOGQ>
bool Z_Q<LOGQ>::is_bigger_than_halfQ()const{
	int length = (LOGQ+63)/64;
	int r = LOGQ -(length-1)*64;
	return (data[length-1]>>(r-1))==1;
}


// convert to big integer of different size
template<int LOGQfrom, int LOGQto>
void resize(const Z_Q<LOGQfrom> &Afrom, Z_Q<LOGQto> &Ato) {
	bool is_negative = Afrom.is_bigger_than_halfQ();
	if(is_negative) {
		Z_Q<LOGQfrom> Afrom_abs(Afrom);
		Afrom_abs.negate();
		Ato.setzero();
		for(int i = 0; i < min(Afrom.get_length(), Ato.get_length()); ++i)
			Ato.data[i] = Afrom_abs.data[i];
		Ato.negate();
	} else {
		Ato.setzero();
		for(int i = 0; i < min(Afrom.get_length(), Ato.get_length()); ++i)
			Ato.data[i] = Afrom.data[i];
	}
	Ato.remove_clutter();
}

template<int LOGQ>
Z_Q<LOGQ>::operator double() const {
	Z_Q<LOGQ> val(*this);
	bool is_negative =  val.is_bigger_than_halfQ();
	if(is_negative)
	    val.negate();
	double result = (double) val[0];
	if(is_negative)
	    result *= -1;
	return result;
}

template<int LOGQ>
int Z_Q<LOGQ>::max_valid_digit() const {
	uint64_t invalid_digit = is_bigger_than_halfQ() ? 0xFFFFFFFFULL : 0;
	for(int i = (LOGQ+63)/64 - 1; i >= 0; --i) {
		if(data[i] != invalid_digit)
			return i;
	}
	return -1;
}