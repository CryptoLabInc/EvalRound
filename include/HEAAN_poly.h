#pragma once

#include "HEAAN.h"

template< int LOGQ, int N, int LOGDELTA >
void eval_poly_deg3( const double c[4],
					 const R_Q_square<   LOGQ         ,N>&   ct,
					 const R_Q_square<2* LOGQ         ,N>&  evk1,
					 const R_Q_square<2*(LOGQ-LOGDELTA),N>& evk2,
						   R_Q_square<   LOGQ-2*LOGDELTA,N>&   p ){
	// bottom level
	R_Q<LOGQ,N> L0[4];
	for(int n=0;n<4;n++){
		L0[n].setzero();
		L0[n][0][0]=uint64_t(c[n]*(1ULL<<LOGDELTA));
		if(n%2==0)
			L0[n][0]*=1ULL<<LOGDELTA;
	}
	// level one
	const R_Q_square<LOGQ,N>& T1=ct;
		  R_Q_square<LOGQ-LOGDELTA,N> L1[2];
	
	for(int n=0; n<2; n++){
		R_Q_square<LOGQ,N> temp; temp=T1; 
		temp    *=L0[2*n+1];
		temp[0] +=L0[2*n  ];
		RS<LOGQ,LOGQ-LOGDELTA,N>(temp,L1[n]);
	}
	//level two
	R_Q_square<LOGQ-LOGDELTA,N> T2;
	{ R_Q_square<LOGQ,N> temp; Mul<LOGQ,N>(T1,T1,evk1,temp); temp*=2;
	  Z_Q<LOGQ> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
	  temp[0][0]-=one;
	  RS<LOGQ, LOGQ-LOGDELTA,N>(temp,T2);
	}
	{R_Q_square<LOGQ-LOGDELTA,N> temp;
	 Mul<LOGQ-LOGDELTA,N>(T2, L1[1], evk2, temp);
	 L1[0]*=1ULL<<LOGDELTA; temp+=L1[0];
	 RS<LOGQ-LOGDELTA, LOGQ-LOGDELTA*2, N>(temp,p);
	}
}

template< int LOGQ, int N, int LOGDELTA >
void eval_poly_deg7( const double c[8],
					 const R_Q_square<   LOGQ            ,N>&   ct,
					 const R_Q_square<2* LOGQ            ,N>&  evk1,
					 const R_Q_square<2*(LOGQ-  LOGDELTA),N>& evk2,
					 const R_Q_square<2*(LOGQ-2*LOGDELTA),N>& evk3,
						   R_Q_square<   LOGQ-3*LOGDELTA ,N>&   p ){
	// bottom level
	R_Q<LOGQ,N> L0[8];
	for(int n=0;n<8;n++){
		L0[n].setzero();
		L0[n][0][0]=uint64_t(c[n]*(1ULL<<LOGDELTA));
		if(n%2==0)
			L0[n][0]*=1ULL<<LOGDELTA;
	}
	// level one
	const R_Q_square<LOGQ,N>& T1=ct;
		  R_Q_square<LOGQ-LOGDELTA  ,N> L1[4];

	for(int n=0; n<4; n++){
		R_Q_square<LOGQ,N> temp; temp=T1; 
		temp    *=L0[2*n+1];
		temp[0] +=L0[2*n  ];
		RS<LOGQ,LOGQ-LOGDELTA,N>(temp,L1[n]);
	}
	//level two
	R_Q_square<LOGQ-LOGDELTA,N> T2;
	{ R_Q_square<LOGQ,N> temp; Mul<LOGQ,N>(T1,T1,evk1,temp); temp*=2;
	  Z_Q<LOGQ> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
	  temp[0][0]-=one;
	  RS<LOGQ, LOGQ-LOGDELTA,N>(temp,T2);
	}
	R_Q_square<LOGQ-LOGDELTA*2,N> L2[2];
	
	for(int n=0; n<2; n++){
		R_Q_square<LOGQ-LOGDELTA ,N> temp;
		Mul<LOGQ-LOGDELTA,N>(T2, L1[2*n+1], evk2, temp);
		L1[2*n]*=1ULL<<LOGDELTA; temp+=L1[2*n];
		RS<LOGQ-LOGDELTA,LOGQ-LOGDELTA*2,N>(temp,L2[n]);
	}

	//level three
	R_Q_square<LOGQ-LOGDELTA*2,N> T4;
	{ 
		R_Q_square<LOGQ-LOGDELTA,N> temp; Mul<LOGQ-LOGDELTA,N>(T2,T2,evk2,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA, LOGQ-LOGDELTA*2,N>(temp,T4);
	}
	T4.print();
	L2[0].print();
	L2[1].print();
	{
		R_Q_square<LOGQ-LOGDELTA*2,N> temp;
		Mul<LOGQ-LOGDELTA*2,N>(T4, L2[1], evk3, temp);
		L2[0]*=1ULL<<LOGDELTA; temp+=L2[0];
		
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3, N>(temp,p);
		
	}
}

template< int LOGQ, int N, int LOGDELTA >
void eval_poly_deg15( const double c[16],
					 const R_Q_square<   LOGQ            ,N>&   ct,
					 const R_Q_square<2* LOGQ            ,N>&  evk1,
					 const R_Q_square<2*(LOGQ-  LOGDELTA),N>& evk2,
					 const R_Q_square<2*(LOGQ-2*LOGDELTA),N>& evk3,
				     const R_Q_square<2*(LOGQ-3*LOGDELTA),N>& evk4,
						   R_Q_square<   LOGQ-4*LOGDELTA ,N>&   p ){
	// bottom level
	R_Q<LOGQ,N> L0[16];
	for(int n=0;n<16;n++){
		L0[n].setzero();
		L0[n][0][0]=uint64_t(c[n]*(1ULL<<LOGDELTA));
		if(n%2==0)
			L0[n][0]*=1ULL<<LOGDELTA;
	}
	// level one
	const R_Q_square<LOGQ,N>& T1=ct;
		  R_Q_square<LOGQ-LOGDELTA  ,N> L1[8];

	for(int n=0; n<8; n++){
		R_Q_square<LOGQ,N> temp; temp=T1; 
		temp    *=L0[2*n+1];
		temp[0] +=L0[2*n  ];
		RS<LOGQ,LOGQ-LOGDELTA,N>(temp,L1[n]);
	}
	//level two
	R_Q_square<LOGQ-LOGDELTA,N> T2;
	{ R_Q_square<LOGQ,N> temp; Mul<LOGQ,N>(T1,T1,evk1,temp); temp*=2;
	  Z_Q<LOGQ> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
	  temp[0][0]-=one;
	  RS<LOGQ, LOGQ-LOGDELTA,N>(temp,T2);
	}
	R_Q_square<LOGQ-LOGDELTA*2,N> L2[4];
	
	for(int n=0; n<4; n++){
		R_Q_square<LOGQ-LOGDELTA ,N> temp;
		Mul<LOGQ-LOGDELTA,N>(T2, L1[2*n+1], evk2, temp);
		L1[2*n]*=1ULL<<LOGDELTA; temp+=L1[2*n];
		RS<LOGQ-LOGDELTA,LOGQ-LOGDELTA*2,N>(temp,L2[n]);
	}

	//level three
	R_Q_square<LOGQ-LOGDELTA*2,N> T4;
	{ 
		R_Q_square<LOGQ-LOGDELTA,N> temp; Mul<LOGQ-LOGDELTA,N>(T2,T2,evk2,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA, LOGQ-LOGDELTA*2,N>(temp,T4);
	}
	R_Q_square<LOGQ-LOGDELTA*3,N> L3[2];

	for(int n=0; n<2; n++){
		R_Q_square<LOGQ-LOGDELTA*2,N> temp;
		Mul<LOGQ-LOGDELTA*2,N>(T4, L2[2*n+1], evk3, temp);
		L2[2*n]*=1ULL<<LOGDELTA; temp+=L2[2*n];
		
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3, N>(temp,L3[n]);
		
	}

	//level four
	R_Q_square<LOGQ-LOGDELTA*3,N> T8;
	{ 
		R_Q_square<LOGQ-LOGDELTA*2,N> temp; Mul<LOGQ-LOGDELTA*2,N>(T4,T4,evk3,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*2> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3,N>(temp,T8);
	}
	
	{
		R_Q_square<LOGQ-LOGDELTA*3,N> temp;
		Mul<LOGQ-LOGDELTA*3,N>(T8, L3[1], evk4, temp);
		L3[0]*=1ULL<<LOGDELTA; temp+=L3[0];
		
		RS<LOGQ-LOGDELTA*3, LOGQ-LOGDELTA*4, N>(temp,p);
		
	}
}


template< int LOGQ, int N, int LOGDELTA >
void eval_poly_deg31( const double c[32],
					 const R_Q_square<   LOGQ            ,N>&   ct,
					 const R_Q_square<2* LOGQ            ,N>&  evk1,
					 const R_Q_square<2*(LOGQ-  LOGDELTA),N>& evk2,
					 const R_Q_square<2*(LOGQ-2*LOGDELTA),N>& evk3,
				     const R_Q_square<2*(LOGQ-3*LOGDELTA),N>& evk4,
					 const R_Q_square<2*(LOGQ-4*LOGDELTA),N>& evk5,
						   R_Q_square<   LOGQ-5*LOGDELTA ,N>&   p ){
	// bottom level
	R_Q<LOGQ,N> L0[32];
	for(int n=0;n<32;n++){
		L0[n].setzero();
		L0[n][0][0]=uint64_t(fabs(c[n])*(1ULL<<LOGDELTA));
		if(n%2==0)
			L0[n][0]*=1ULL<<LOGDELTA;
		if(c[n]<0) L0[n][0].negate();
	}
	L0[0].print_unsigned();
	L0[1].print_unsigned();
	L0[2].print_unsigned();
	L0[3].print_unsigned();
	// level one
	const R_Q_square<LOGQ,N>& T1=ct;
		  R_Q_square<LOGQ-LOGDELTA  ,N> L1[16];
	T1.print();
	for(int n=0; n<16; n++){
		R_Q_square<LOGQ,N> temp; temp=T1; 
		temp    *=L0[2*n+1];
		temp[0] +=L0[2*n  ];
		RS<LOGQ,LOGQ-LOGDELTA,N>(temp,L1[n]);
	}
	L1[0].print();
	L1[1].print();
	//level two
	R_Q_square<LOGQ-LOGDELTA,N> T2;
	{ R_Q_square<LOGQ,N> temp; Mul<LOGQ,N>(T1,T1,evk1,temp); temp*=2;
	  Z_Q<LOGQ> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
	  temp[0][0]-=one;
	  RS<LOGQ, LOGQ-LOGDELTA,N>(temp,T2);
	}
	T2.print();
	R_Q_square<LOGQ-LOGDELTA*2,N> L2[8];
	
	for(int n=0; n<8; n++){
		R_Q_square<LOGQ-LOGDELTA ,N> temp;
		Mul<LOGQ-LOGDELTA,N>(T2, L1[2*n+1], evk2, temp);
		L1[2*n]*=1ULL<<LOGDELTA; temp+=L1[2*n];
		RS<LOGQ-LOGDELTA,LOGQ-LOGDELTA*2,N>(temp,L2[n]);
	}
	L2[0].print();
	//level three
	R_Q_square<LOGQ-LOGDELTA*2,N> T4;
	{ 
		R_Q_square<LOGQ-LOGDELTA,N> temp; Mul<LOGQ-LOGDELTA,N>(T2,T2,evk2,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA, LOGQ-LOGDELTA*2,N>(temp,T4);
	}
	R_Q_square<LOGQ-LOGDELTA*3,N> L3[4];

	for(int n=0; n<4; n++){
		R_Q_square<LOGQ-LOGDELTA*2,N> temp;
		Mul<LOGQ-LOGDELTA*2,N>(T4, L2[2*n+1], evk3, temp);
		L2[2*n]*=1ULL<<LOGDELTA; temp+=L2[2*n];
		
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3, N>(temp,L3[n]);
		
	}
	L3[0].print();
	//level four
	R_Q_square<LOGQ-LOGDELTA*3,N> T8;
	{ 
		R_Q_square<LOGQ-LOGDELTA*2,N> temp; Mul<LOGQ-LOGDELTA*2,N>(T4,T4,evk3,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*2> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3,N>(temp,T8);
	}
	R_Q_square<LOGQ-LOGDELTA*4,N> L4[2];

	for(int n=0; n<2; n++){
		R_Q_square<LOGQ-LOGDELTA*3,N> temp;
		Mul<LOGQ-LOGDELTA*3,N>(T8, L3[2*n+1], evk4, temp);
		L3[2*n]*=1ULL<<LOGDELTA; temp+=L3[2*n];
		
		RS<LOGQ-LOGDELTA*3, LOGQ-LOGDELTA*4, N>(temp,L4[n]);
		
	}
	L4[0].print();
	
	//level five
	R_Q_square<LOGQ-LOGDELTA*4,N> T16;
	{ 
		R_Q_square<LOGQ-LOGDELTA*3,N> temp; Mul<LOGQ-LOGDELTA*3,N>(T8,T8,evk4,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*3> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*3, LOGQ-LOGDELTA*4,N>(temp,T16);
	}
	
	{
		R_Q_square<LOGQ-LOGDELTA*4,N> temp;
		Mul<LOGQ-LOGDELTA*4,N>(T16, L4[1], evk5, temp);
		L4[0]*=1ULL<<LOGDELTA; temp+=L4[0];
		
		RS<LOGQ-LOGDELTA*4, LOGQ-LOGDELTA*5, N>(temp,p);
		
	}
}


template< int LOGQ, int N, int LOGDELTA >
void eval_poly_deg63( const double c[64],
					 const R_Q_square<   LOGQ            ,N>&   ct,
					 const R_Q_square<2* LOGQ            ,N>&  evk1,
					 const R_Q_square<2*(LOGQ-  LOGDELTA),N>& evk2,
					 const R_Q_square<2*(LOGQ-2*LOGDELTA),N>& evk3,
				     const R_Q_square<2*(LOGQ-3*LOGDELTA),N>& evk4,
					 const R_Q_square<2*(LOGQ-4*LOGDELTA),N>& evk5,
					 const R_Q_square<2*(LOGQ-5*LOGDELTA),N>& evk6,
						   R_Q_square<   LOGQ-6*LOGDELTA ,N>&   p ){
	// bottom level
	R_Q<LOGQ,N> L0[64];
	for(int n=0;n<64;n++){
		L0[n].setzero();
		L0[n][0][0]=uint64_t(c[n]*(1ULL<<LOGDELTA));
		if(n%2==0)
			L0[n][0]*=1ULL<<LOGDELTA;
	}
	// level one
	const R_Q_square<LOGQ,N>& T1=ct;
		  R_Q_square<LOGQ-LOGDELTA  ,N> L1[32];

	for(int n=0; n<32; n++){
		R_Q_square<LOGQ,N> temp; temp=T1; 
		temp    *=L0[2*n+1];
		temp[0] +=L0[2*n  ];
		RS<LOGQ,LOGQ-LOGDELTA,N>(temp,L1[n]);
	}
	//level two
	R_Q_square<LOGQ-LOGDELTA,N> T2;
	{ R_Q_square<LOGQ,N> temp; Mul<LOGQ,N>(T1,T1,evk1,temp); temp*=2;
	  Z_Q<LOGQ> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
	  temp[0][0]-=one;
	  RS<LOGQ, LOGQ-LOGDELTA,N>(temp,T2);
	}
	R_Q_square<LOGQ-LOGDELTA*2,N> L2[16];
	
	for(int n=0; n<16; n++){
		R_Q_square<LOGQ-LOGDELTA ,N> temp;
		Mul<LOGQ-LOGDELTA,N>(T2, L1[2*n+1], evk2, temp);
		L1[2*n]*=1ULL<<LOGDELTA; temp+=L1[2*n];
		RS<LOGQ-LOGDELTA,LOGQ-LOGDELTA*2,N>(temp,L2[n]);
	}

	//level three
	R_Q_square<LOGQ-LOGDELTA*2,N> T4;
	{ 
		R_Q_square<LOGQ-LOGDELTA,N> temp; Mul<LOGQ-LOGDELTA,N>(T2,T2,evk2,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA, LOGQ-LOGDELTA*2,N>(temp,T4);
	}
	R_Q_square<LOGQ-LOGDELTA*3,N> L3[8];

	for(int n=0; n<8; n++){
		R_Q_square<LOGQ-LOGDELTA*2,N> temp;
		Mul<LOGQ-LOGDELTA*2,N>(T4, L2[2*n+1], evk3, temp);
		L2[2*n]*=1ULL<<LOGDELTA; temp+=L2[2*n];
		
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3, N>(temp,L3[n]);
		
	}

	//level four
	R_Q_square<LOGQ-LOGDELTA*3,N> T8;
	{ 
		R_Q_square<LOGQ-LOGDELTA*2,N> temp; Mul<LOGQ-LOGDELTA*2,N>(T4,T4,evk3,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*2> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3,N>(temp,T8);
	}
	R_Q_square<LOGQ-LOGDELTA*4,N> L4[4];

	for(int n=0; n<4; n++){
		R_Q_square<LOGQ-LOGDELTA*3,N> temp;
		Mul<LOGQ-LOGDELTA*3,N>(T8, L3[2*n+1], evk4, temp);
		L3[2*n]*=1ULL<<LOGDELTA; temp+=L3[2*n];
		
		RS<LOGQ-LOGDELTA*3, LOGQ-LOGDELTA*4, N>(temp,L4[n]);
		
	}

	
	//level five
	R_Q_square<LOGQ-LOGDELTA*4,N> T16;
	{ 
		R_Q_square<LOGQ-LOGDELTA*3,N> temp; Mul<LOGQ-LOGDELTA*3,N>(T8,T8,evk4,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*3> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*3, LOGQ-LOGDELTA*4,N>(temp,T16);
	}
	R_Q_square<LOGQ-LOGDELTA*5,N> L5[2];

	for(int n=0; n<2; n++){
		R_Q_square<LOGQ-LOGDELTA*4,N> temp;
		Mul<LOGQ-LOGDELTA*4,N>(T16, L4[2*n+1], evk5, temp);
		L4[2*n]*=1ULL<<LOGDELTA; temp+=L4[2*n];
		
		RS<LOGQ-LOGDELTA*4, LOGQ-LOGDELTA*5, N>(temp,L5[n]);
		
	}

	//level six
	R_Q_square<LOGQ-LOGDELTA*5,N> T32;
	{ 
		R_Q_square<LOGQ-LOGDELTA*4,N> temp; Mul<LOGQ-LOGDELTA*4,N>(T16,T16,evk5,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*4> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*4, LOGQ-LOGDELTA*5,N>(temp,T32);
	}
	
	{	R_Q_square<LOGQ-LOGDELTA*5,N> temp;
		Mul<LOGQ-LOGDELTA*5,N>(T32, L5[1], evk6, temp);
		L5[0]*=1ULL<<LOGDELTA; temp+=L5[0];
		
		RS<LOGQ-LOGDELTA*5, LOGQ-LOGDELTA*6, N>(temp,p);
		
	}
}

template< int LOGQ, int N, int LOGDELTA >
void eval_poly_deg127( const double c[128],
					 const R_Q_square<   LOGQ            ,N>&   ct,
					 const R_Q_square<2* LOGQ            ,N>&  evk1,
					 const R_Q_square<2*(LOGQ-  LOGDELTA),N>& evk2,
					 const R_Q_square<2*(LOGQ-2*LOGDELTA),N>& evk3,
				     const R_Q_square<2*(LOGQ-3*LOGDELTA),N>& evk4,
					 const R_Q_square<2*(LOGQ-4*LOGDELTA),N>& evk5,
					 const R_Q_square<2*(LOGQ-5*LOGDELTA),N>& evk6,
					 const R_Q_square<2*(LOGQ-6*LOGDELTA),N>& evk7,
						   R_Q_square<   LOGQ-7*LOGDELTA ,N>&   p ){
	// bottom level
	R_Q<LOGQ,N> L0[128];
	for(int n=0;n<128;n++){
		L0[n].setzero();
		L0[n][0][0]=uint64_t(c[n]*(1ULL<<LOGDELTA));
		if(n%2==0)
			L0[n][0]*=1ULL<<LOGDELTA;
	}
	// level one
	const R_Q_square<LOGQ,N>& T1=ct;
		  R_Q_square<LOGQ-LOGDELTA  ,N> L1[64];

	for(int n=0; n<64; n++){
		R_Q_square<LOGQ,N> temp; temp=T1; 
		temp    *=L0[2*n+1];
		temp[0] +=L0[2*n  ];
		RS<LOGQ,LOGQ-LOGDELTA,N>(temp,L1[n]);
	}
	//level two
	R_Q_square<LOGQ-LOGDELTA,N> T2;
	{ R_Q_square<LOGQ,N> temp; Mul<LOGQ,N>(T1,T1,evk1,temp); temp*=2;
	  Z_Q<LOGQ> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
	  temp[0][0]-=one;
	  RS<LOGQ, LOGQ-LOGDELTA,N>(temp,T2);
	}
	R_Q_square<LOGQ-LOGDELTA*2,N> L2[32];
	
	for(int n=0; n<32; n++){
		R_Q_square<LOGQ-LOGDELTA ,N> temp;
		Mul<LOGQ-LOGDELTA,N>(T2, L1[2*n+1], evk2, temp);
		L1[2*n]*=1ULL<<LOGDELTA; temp+=L1[2*n];
		RS<LOGQ-LOGDELTA,LOGQ-LOGDELTA*2,N>(temp,L2[n]);
	}

	//level three
	R_Q_square<LOGQ-LOGDELTA*2,N> T4;
	{ 
		R_Q_square<LOGQ-LOGDELTA,N> temp; Mul<LOGQ-LOGDELTA,N>(T2,T2,evk2,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA, LOGQ-LOGDELTA*2,N>(temp,T4);
	}
	R_Q_square<LOGQ-LOGDELTA*3,N> L3[16];

	for(int n=0; n<16; n++){
		R_Q_square<LOGQ-LOGDELTA*2,N> temp;
		Mul<LOGQ-LOGDELTA*2,N>(T4, L2[2*n+1], evk3, temp);
		L2[2*n]*=1ULL<<LOGDELTA; temp+=L2[2*n];
		
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3, N>(temp,L3[n]);
		
	}

	//level four
	R_Q_square<LOGQ-LOGDELTA*3,N> T8;
	{ 
		R_Q_square<LOGQ-LOGDELTA*2,N> temp; Mul<LOGQ-LOGDELTA*2,N>(T4,T4,evk3,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*2> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*2, LOGQ-LOGDELTA*3,N>(temp,T8);
	}
	R_Q_square<LOGQ-LOGDELTA*4,N> L4[8];

	for(int n=0; n<8; n++){
		R_Q_square<LOGQ-LOGDELTA*3,N> temp;
		Mul<LOGQ-LOGDELTA*3,N>(T8, L3[2*n+1], evk4, temp);
		L3[2*n]*=1ULL<<LOGDELTA; temp+=L3[2*n];
		
		RS<LOGQ-LOGDELTA*3, LOGQ-LOGDELTA*4, N>(temp,L4[n]);
		
	}

	
	//level five
	R_Q_square<LOGQ-LOGDELTA*4,N> T16;
	{ 
		R_Q_square<LOGQ-LOGDELTA*3,N> temp; Mul<LOGQ-LOGDELTA*3,N>(T8,T8,evk4,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*3> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*3, LOGQ-LOGDELTA*4,N>(temp,T16);
	}
	R_Q_square<LOGQ-LOGDELTA*5,N> L5[4];

	for(int n=0; n<4; n++){
		R_Q_square<LOGQ-LOGDELTA*4,N> temp;
		Mul<LOGQ-LOGDELTA*4,N>(T16, L4[2*n+1], evk5, temp);
		L4[2*n]*=1ULL<<LOGDELTA; temp+=L4[2*n];
		
		RS<LOGQ-LOGDELTA*4, LOGQ-LOGDELTA*5, N>(temp,L5[n]);
		
	}

	//level six
	R_Q_square<LOGQ-LOGDELTA*5,N> T32;
	{ 
		R_Q_square<LOGQ-LOGDELTA*4,N> temp; Mul<LOGQ-LOGDELTA*4,N>(T16,T16,evk5,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*4> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*4, LOGQ-LOGDELTA*5,N>(temp,T32);
	}
	R_Q_square<LOGQ-LOGDELTA*6,N> L6[2];

	for(int n=0; n<2; n++){
		R_Q_square<LOGQ-LOGDELTA*5,N> temp;
		Mul<LOGQ-LOGDELTA*5,N>(T32, L5[2*n+1], evk6, temp);
		L5[2*n]*=1ULL<<LOGDELTA; temp+=L5[2*n];
		
		RS<LOGQ-LOGDELTA*5, LOGQ-LOGDELTA*6, N>(temp,L6[n]);
		
	}

	//level seven
	R_Q_square<LOGQ-LOGDELTA*6,N> T64;
	{ 
		R_Q_square<LOGQ-LOGDELTA*5,N> temp; Mul<LOGQ-LOGDELTA*5,N>(T32,T32,evk6,temp); temp*=2;
		Z_Q<LOGQ-LOGDELTA*5> one; one.setzero(); one[0]=1ULL<<LOGDELTA; one*=1ULL<<LOGDELTA; 
		temp[0][0]-=one;
		RS<LOGQ-LOGDELTA*5, LOGQ-LOGDELTA*6,N>(temp,T64);
	}
	
	{	R_Q_square<LOGQ-LOGDELTA*6,N> temp;
		Mul<LOGQ-LOGDELTA*6,N>(T64, L6[1], evk7, temp);
		L6[0]*=1ULL<<LOGDELTA; temp+=L6[0];
		
		RS<LOGQ-LOGDELTA*6, LOGQ-LOGDELTA*7, N>(temp,p);
		
	}
}


void convert_frc_tou( const double* c, double* u, const int N ){
	if(N==4){ u[0] = c[0];
			  u[1] = c[1]-c[3];
			  u[2] = c[2];
			  u[3] = c[3]*2;}
	else    { double* c0=new double[N/2]; for(int i=0; i<N/2; i++) c0[i] = c[i    ];
			  double* c1=new double[N/2]; for(int i=0; i<N/2; i++) c1[i] = c[i+N/2];
			  for(int i=1; i<N/2; i++){
				  c0[N/2-i]-=c1[i]; c1[i]*=2;}
			  convert_frc_tou(c0,u    ,N/2);
			  convert_frc_tou(c1,u+N/2,N/2);
			  delete[] c0;
			  delete[] c1;
	}
}