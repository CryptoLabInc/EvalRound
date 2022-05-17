#pragma once

#include <math.h>
#include "HEAAN_matrix.h"
#include "HEAAN_poly.h"


template<int LOGQ, int LOGN>
void CoeffToSlot( const R_Q_square<  LOGQ,1<<LOGN>& ct,
				  const R_Q_square<2*LOGQ,1<<LOGN>& rkey,
				  const R_Q_square<2*LOGQ,1<<LOGN>& ckey,
						R_Q_square<  LOGQ,1<<LOGN> ct_[2]){
	const int N = 1 << LOGN;
	double U0r[N/2][N/2];
	double U0i[N/2][N/2]; get_U0<N>(U0r, U0i);
	double Ar[N/2][N/2];
	double Ai[N/2][N/2];

	for(int i=0; i<N/2; i++)
	for(int j=0; j<N/2; j++){
		Ar[i][j] =  U0r[j][i]/N;
		Ai[i][j] = -U0i[j][i]/N;
	}
	R_Q_square<LOGQ,N> ct1; linear_transform<LOGQ,LOGN,50>(Ar,Ai,ct,rkey,ct1);
	R_Q_square<LOGQ,N> ct2; conj(ct1,ckey,ct2);

	R_Q<LOGQ,N> pti; pti.setzero();
	pti[N/2][0] = 1; pti[N/2].negate();
	ct_[1] = ct1; ct_[1] -=ct2; ct_[1]*=pti;
	ct_[0] = ct1; ct_[0] +=ct2;
}

template<int LOGQ, int LOGN, int LOGDELTA, int G = 1>
void CoeffToSlot (	const R_Q_square<  LOGQ,1<<LOGN>& ct,
					const int s[1 << LOGN],
					R_Q_square<  LOGQ,1<<LOGN> ct_[2]){
	const int N = 1 << LOGN;
	static bool is_init = false;
	static SparseDiagonal<N/2,3> U0r[LOGN-1];
    static SparseDiagonal<N/2,3> U0i[LOGN-1];

	if(!is_init) {
		splitU0NR<LOGN>(U0r, U0i);
    	for (int n = 0; n < LOGN - 1; n++) {
			U0r[n].transpose();
			U0i[n].transpose();
			U0i[n].negate();
			U0r[n] *= (n == 0) ? 0.25 : 0.5;
			U0i[n] *= (n == 0) ? 0.25 : 0.5;
		}
		is_init = true;
	}

	R_Q_square<LOGQ, N> U0ct;
	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA, 3, LOGN-1, G>(U0r, U0i, ct, s, U0ct);

	R_Q_square<LOGQ, N> U0ct_conj;
	int s_conj[N];
	R_Q_square<2*LOGQ, 1 << LOGN> ckey;
	conj<N>(s, s_conj);
	HEAAN<LOGQ,N>::swkgen(s_conj, s, ckey);
	conj(U0ct, ckey, U0ct_conj);

	R_Q<LOGQ,N> pti; pti.setzero();
	pti[N/2][0] = 1; pti[N/2].negate();
	ct_[0] = U0ct; ct_[0] += U0ct_conj;
	ct_[1] = U0ct; ct_[1] -= U0ct_conj; ct_[1]*=pti;
}

//---------------------------------------------
// CoeffToSlot with splitU0 with N=16
//---------------------------------------------
template<int LOGQ>
void CoeffToSlot(const R_Q_square<  LOGQ, 16>& ct,
				 const R_Q_square<2 * LOGQ, 16>& rkey1,
				 const R_Q_square<2 * LOGQ, 16>& rkey2,
				 const R_Q_square<2 * LOGQ, 16>& rkey4,
				 const R_Q_square<2 * LOGQ, 16>& rkey6,
				 const R_Q_square<2 * LOGQ, 16>& rkey7,
				 const R_Q_square<2 * LOGQ, 16>& ckey,
					   R_Q_square<  LOGQ, 16> ct_[2]) {
	SparseDiagonal<8, 3> U0r[3];
	SparseDiagonal<8, 3> U0i[3];
	splitU0NR<4>(U0r, U0i);

	for (int n = 0; n < 3; n++) {
		U0r[n].transpose();
		U0i[n].transpose();
		U0i[n].negate();
	}
	U0r[0] *= (1. / 16);
	U0i[0] *= (1. / 16);

	const R_Q_square<2 * LOGQ, 16>* rkey[3];
	const R_Q_square<    LOGQ, 16> ct1, ct2, temp;
	rkey[0] = &rkey1;
	rkey[1] = &rkey1;
	rkey[2] = &rkey1;
	linear_transform<LOGQ, 16, 50, 3>(U0r[0], U0i[0], ct, rkey[3], ct1);

	rkey[0] = &rkey2;
	rkey[1] = &rkey2;
	rkey[2] = &rkey6;
	linear_transform<LOGQ, 16, 50, 3>(U0r[1], U0i[1], ct1, rkey[3], temp);
		
	rkey[0] = &rkey1;
	rkey[1] = &rkey1;
	rkey[2] = &rkey7;
	linear_transform<LOGQ, 16, 50, 3>(U0r[2], U0i[2], temp, rkey[3], ct1);

	conj<LOGQ, 16>(ct1, ckey, ct2);
	R_Q<LOGQ, 16>pti; pti.setzero();
	pti[8][0] = 1; pti[8].negate();
	ct_[1] = ct1; ct_[1] -= ct2;
	ct_[0] = ct1; ct_[0] -= ct2;

}


template<int LOGQ, int LOGN>
void SlotToCoeff( const R_Q_square<  LOGQ,1<<LOGN>&  ct0,
				  const R_Q_square<  LOGQ,1<<LOGN>&  ct1,
				  const R_Q_square<2*LOGQ,1<<LOGN>& rkey,
					    R_Q_square<  LOGQ,1<<LOGN>& ct_ ){
	const int N = 1 << LOGN;
	double U0r[N/2][N/2];
	double U0i[N/2][N/2]; get_U0<N>(U0r,U0i);
	linear_transform<LOGQ,LOGN,50>(U0r,U0i,ct0,rkey,ct_);

	double iU0r[N/2][N/2];
	double iU0i[N/2][N/2];

	for(int i=0; i<N/2; i++)
	for(int j=0; j<N/2; j++){
		iU0r[i][j] =-U0i[i][j];
		iU0i[i][j] = U0r[i][j];
	}
	R_Q_square<LOGQ,N> temp;
	linear_transform<LOGQ,LOGN,50>(iU0r,iU0i,ct1,rkey,temp);
	ct_+=temp;
}

template<int LOGQ, int LOGN, int LOGDELTA, int G>
void SlotToCoeff( const R_Q_square<  LOGQ,1<<LOGN>&  ct0,
				  const R_Q_square<  LOGQ,1<<LOGN>&  ct1,
				  const int skey[1<<LOGN],
					    R_Q_square<  LOGQ,1<<LOGN>& ct_ ){	
	const int N = 1 << LOGN;
	static bool is_init = false;
	static SparseDiagonal<N/2,3> U0r[LOGN-1];
    static SparseDiagonal<N/2,3> U0i[LOGN-1];
	static SparseDiagonal<N/2,3> iU0r[LOGN-1];
    static SparseDiagonal<N/2,3> iU0i[LOGN-1];

	if(!is_init) {
		SparseDiagonal<N/2,3> U0r_temp[LOGN-1];
    	SparseDiagonal<N/2,3> U0i_temp[LOGN-1];

		splitU0NR<LOGN>(U0r_temp, U0i_temp);
		for(int i = 0; i < LOGN - 1; ++i) {
			U0r[i] = U0r_temp[LOGN - 2 - i];
			U0i[i] = U0i_temp[LOGN - 2 - i];
		}

		for(int i = 0; i < LOGN - 1; ++i) {
			if(i == 0) {
				iU0r[i] = U0i_temp[LOGN - 2 - i];
				iU0r[i].negate();
				iU0i[i] = U0r_temp[LOGN - 2 - i];
			} else {
				iU0r[i] = U0r_temp[LOGN - 2 - i];
				iU0i[i] = U0i_temp[LOGN - 2 - i];
			}
		}		
		is_init = true;
	}

	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA, 3, LOGN-1, G>(U0r, U0i, ct0, skey, ct_);
	R_Q_square<LOGQ, N> ct2;
	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA, 3, LOGN-1, G>(iU0r, iU0i, ct1, skey, ct2);
	ct_ += ct2;
}


template<int LOGN>
void split_U0( double Ar[LOGN][1<<(LOGN-1)][1<<(LOGN-1)],
			   double Ai[LOGN][1<<(LOGN-1)][1<<(LOGN-1)] ){
	int N = 1<<LOGN;
	for(int n=0; n<LOGN; n++)
	for(int i=0; i<N/2;  i++)
	for(int j=0; j<N/2;  j++){ Ar[n][i][j] = 0;
							   Ai[n][i][j] = 0;}
	for(int n=0; n<LOGN-1; n++){
		int M = N>>(n+2);
		for(int i=0, fivei=1<<n; i<M; i++, fivei = (fivei*5)%(2*N)){
			Ar[n][i][i  ] = 1; Ar[n][i+M][i] = 1;
			Ar[n][i][i+M] = cos(PI/N*fivei); Ar[n][i+M][i+M] = -Ar[n][i][i+M];
			Ai[n][i][i+M] = sin(PI/N*fivei); Ai[n][i+M][i+M] = -Ai[n][i][i+M];
		}									   
		for(int k=1; k<(1<<n); k++){ int off = k*2*M;
		for(int i=0; i<M; i++){
			Ar[n][i  +off][i  +off] = Ar[n][i  ][i  ];
			Ar[n][i  +off][i+M+off] = Ar[n][i  ][i+M];
			Ar[n][i+M+off][i  +off] = Ar[n][i+M][i  ];
			Ar[n][i+M+off][i+M+off] = Ar[n][i+M][i+M];
			Ai[n][i  +off][i+M+off] = Ai[n][i  ][i+M];
			Ai[n][i+M+off][i+M+off] = Ai[n][i+M][i+M];
		}}
	}
	for(int i=0; i<N/2; i++){
		int j=0; int i_=i;
		for(int k=0; k<LOGN; k++){
			j+= (i_&1)<<(LOGN-2-k); 
			i_>>=1;
		}

		Ar[LOGN-1][i][j] = 1;
	}
}

template<int LOGQ, int N, int LOGDELTA, int K>
void EvalMod( const R_Q_square<LOGQ, N>& ct,
		      const R_Q_square<2*(LOGQ-   LOGDELTA), N>& evk1,  // poly
		      const R_Q_square<2*(LOGQ- 2*LOGDELTA), N>& evk2,  // poly
		      const R_Q_square<2*(LOGQ- 3*LOGDELTA), N>& evk3,  // poly
		      const R_Q_square<2*(LOGQ- 4*LOGDELTA), N>& evk4,  // poly
		      const R_Q_square<2*(LOGQ- 5*LOGDELTA), N>& evk5,  // poly
		      const R_Q_square<2*(LOGQ- 6*LOGDELTA), N>& evk6,  // poly
		      const R_Q_square<2*(LOGQ- 7*LOGDELTA), N>& evk7,  // poly
		      const R_Q_square<2*(LOGQ- 8*LOGDELTA), N>& evk8,  // double angle
		      const R_Q_square<2*(LOGQ- 9*LOGDELTA), N>& evk9,  // double angle
		      const R_Q_square<2*(LOGQ-10*LOGDELTA), N>& evk10, // arcsine
		      const R_Q_square<2*(LOGQ-11*LOGDELTA), N>& evk11, // arcsine
			        R_Q_square<  (LOGQ-12*LOGDELTA), N>& p    )
{
   	// ct1 = ct * (1/K)
    R_Q_square<LOGQ,N> ct_copy=ct; ct_copy*=uint64_t((1ULL<<LOGDELTA)/static_cast<double>(K));
	R_Q_square<LOGQ-LOGDELTA,N> ct1;
	RS<LOGQ,LOGQ-LOGDELTA,N>(ct_copy, ct1);

    // shift : ct1 -= (1/(4K))	
	Z_Q<LOGQ-LOGDELTA> qtK; qtK.setzero(); qtK[0]=((1ULL<<LOGDELTA)/static_cast<double>(4 * K));
	ct1[0][0] -= qtK;
	
	// Evaluate cos(12 * pi * t) for t in [-1, 1]
	const double u[128] = { 0.0915790575476516607020301925599,
	0.965641566848845149012510534562e-59,
	-0.912759680952499166357451660022,
    0.72729008060602405199483586801e-60,
	0.683392580248665849727127512198,
	0.154791006261887519671458670466e-58,
	0.575775246994246272241797285215,
	-0.136985563540168874690052022864e-58,
	0.187114138152947106022653825575,
	0.102650106470327625507668303615e-58,
	-0.504297367752770601815257078595,
	0.486551958999121159704498347834e-59,
	0.289005559335286254998214282896,
	0.200070598136902342290684628231e-58,
	0.694196477630006039024795353227,
	-0.349522721172103271129549737324e-58,
	-0.242851623141933344719165890514,
	0.136029230370127877805335898372e-58,
	-2.30916607014755066824710202036,
	0.994380409496633126454818333524e-59,
	1.15407689315779284424683872544,
	0.185641029151842855675889943975e-58,
	1.06065310745523018876903064942,
	-0.406560084721113431774768336533e-58,
	0.15051560693828324212142477633,
	0.205840979436728204601610643516e-58,
	-2.41746339029848750634517853491,
	0.128120176730898416670909833799e-58,
	-2.43751736620331230791270329623,
	-0.114539856917654895989595125497e-58,
	5.15307624783564099246348603548,
	-0.33848810450038771462752191236e-58,
	0.0986801315005327533579026611641,
	0.132729288165145015536869774182e-59,
	-0.345468155933006903000331541412,
	0.13578343133985082423102843549e-58,
	0.729151120285267656632735758534,
	0.250013991159788306629186438875e-59,
	-0.795774692329250253203992479218,
	-0.931219426559138119882490846425e-59,
	0.236142830438603593914802647056,
	-0.203741667517328231145028465474e-60,
	-0.16725951439523069822052223604,
	0.326514751051914207721321646764e-58,
	0.0586275712435874653768699211629,
	-0.45471089433299006342598413292e-59,
	-0.0306752689318116060921053799577,
	-0.142366448958522394340589231042e-58,
	0.00192769568901090841886440702691,
	-0.488121731470298000575689254316e-59,
	-0.000785836600946135926067132990438,
	0.118596480963411712296240204552e-58,
	0.000151780007088732069143467138273,
	-0.343768280944662268680930861256e-59,
	-0.500091787382133142651789122259e-4,
	0.249826728880126386495571834874e-58,
	0.38879047610595341122844730393e-5,
	-0.152666401755828603630105166395e-58,
	-0.106180103005540407481893894173e-5,
	-0.155682523598299607635108926542e-58,
	0.136676998903497224402162705731e-6,
	-0.790352029471194251711516833032e-59,
	-0.315505626933620839817093837299e-7,
	-0.148848126475751357391340497392e-58,
	0.214756106143285000557246040321e-9,
	0.536390856802231900990514395321e-59,
	-0.42548806761554568832727652473e-10,
	0.815114145433551422697688653415e-59,
	0.398324976488989589258445432485e-11,
	0.9611333399876254130393044754e-60,
	-0.685914804932498240770203674461e-12,
	-0.17670731044164017398051957779e-58,
	0.279620038414403592993665155674e-13,
	0.305090657278578534681195546566e-60,
	-0.422890001479951589127860746693e-14,
	0.168645983162959362055132164925e-58,
	0.303447662971201151505863162831e-15,
	0.23653274932686362517006880282e-58,
	-0.406661789712370933678410160068e-16,
	-0.153619696965004903185779892563e-58,
	0.647813085218598664401929532027e-18,
	0.498070252702471503532106762163e-59,
	-0.775198400058309223053130021303e-19,
	0.149844613151333990016655378843e-58,
	0.44196681050341771562411835721e-20,
	-0.623505896275676855402963478954e-59,
	-0.475410307298060468129337639151e-21,
	-0.271267930265432733400273573741e-58,
	0.122064446001408267079280056034e-22,
	0.106909481041688967963405927537e-58,
	-0.118723940204738596543958617601e-23,
	0.892310463489543231557017435125e-59,
	0.552285512058246375349611897523e-25,
	-0.767860086406231809424174035768e-59,
	-0.488263794658988003203379473087e-26,
	-0.186130100084505534524888709035e-58,
	0.258519721546982001245624287151e-28,
	0.62286001111692798936428134889e-59,
	-0.208717637953966088305028656679e-29,
	-0.32824510761745634898787502611e-59,
	0.808675498749480381496271982264e-31,
	0.108537170321480044497171702034e-58,
	-0.598752227930039237595031250652e-32,
	-0.239352521045521854600279092857e-58,
	0.10653898760975395746841119585e-33,
	-0.859435937398968676343777716857e-59,
	-0.726189125961755270810675450126e-35,
	0.159269589359434119035519512998e-58,
	0.238250600694213057857702542552e-36,
	0.278346301852733262725121202673e-58,
	-0.150022779902779568923499691185e-37,
	-0.949645965820337683290068597142e-59,
	0.11382816882153704900481501814e-39,
	0.138105990379003129866545733083e-58,
	-0.66426348943571055232797403343e-41,
	-0.951816762764006819136325330609e-59,
	0.18707199197440726330595479964e-42,
	0.321063835439649063098893913726e-58,
	-0.101470651330168920621632563397e-43,
	-0.388828867362704297152874882769e-58,
	0.132964309476343879769218890712e-45,
	0.262135167646277306036072219519e-58,
	-0.672175309441869954743519075739e-47,
	-0.362324728828612142073378625117e-58,
	0.164364057102755835610123118068e-48,
	-0.113587661494912715799544136381e-58,
	-0.777137584964964435504744978855e-50,
	0.631445225239568037148046872981e-59 }; 
    R_Q_square<LOGQ-8*LOGDELTA,N> ct2;
	eval_poly_deg127<LOGQ-LOGDELTA,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,evk6,evk7,ct2);

	// double angle 1
	R_Q_square<LOGQ-8*LOGDELTA,N> temp1; Mul<LOGQ-8*LOGDELTA,N>(ct2,ct2,evk8,temp1); temp1*=2;
    Z_Q<LOGQ-8*LOGDELTA> one1; one1.setzero(); one1[0]=1ULL<<LOGDELTA; one1*=1ULL<<LOGDELTA;
    temp1[0][0]-=one1;
	R_Q_square<LOGQ-9*LOGDELTA,N> ct3;
	RS<LOGQ-8*LOGDELTA,LOGQ-9*LOGDELTA,N>(temp1, ct3);
	
	// double angle 2
	R_Q_square<LOGQ-9*LOGDELTA,N> temp2; Mul<LOGQ-9*LOGDELTA,N>(ct3,ct3,evk9,temp2); temp2*=2;
    Z_Q<LOGQ-9*LOGDELTA> one2; one2.setzero(); one2[0]=1ULL<<LOGDELTA; one2*=1ULL<<LOGDELTA;
    temp2[0][0]-=one2;
	R_Q_square<LOGQ-10*LOGDELTA,N> ct4;
	RS<LOGQ-9*LOGDELTA,LOGQ-10*LOGDELTA,N>(temp2, ct4);

	// arcsine of degree 3
	// ct5 = x^2 + 6
	R_Q_square<LOGQ-10*LOGDELTA,N> temp3; Mul<LOGQ-10*LOGDELTA,N>(ct4,ct4,evk10,temp3);
	Z_Q<LOGQ-10*LOGDELTA> six; six.setzero(); six[0]=(1ULL<<LOGDELTA)*6; six*=1ULL<<LOGDELTA;
	temp3[0][0]+=six;
	R_Q_square<LOGQ-11*LOGDELTA,N> ct5;
	RS<LOGQ-10*LOGDELTA,LOGQ-11*LOGDELTA,N>(temp3,ct5);
	// ct6 = (1/12pi) * x
	temp3=ct4; temp3 *= static_cast<uint64_t>((1ULL<<LOGDELTA) / (12 * PI));
	R_Q_square<LOGQ-11*LOGDELTA,N> ct6;
	RS<LOGQ-10*LOGDELTA,LOGQ-11*LOGDELTA,N>(temp3,ct6);
	// ct7 = (1/2pi) * (x + (1/6) * x^3)
	R_Q_square<LOGQ-11*LOGDELTA,N> ct7; 
	Mul<LOGQ-11*LOGDELTA,N>(ct5,ct6,evk11,ct7);
	RS<LOGQ-11*LOGDELTA,LOGQ-12*LOGDELTA,N>(ct7,p);	    
}

// Rescale by the size of q should not be performed at the beginning of EvalMod_Kx.
/*
template<int LOGQ, int N, int LOGq, int LOGDELTA>
void EvalMod_K4( const R_Q_square<  LOGQ,N>& ct,
				  const R_Q_square<2*(LOGQ-LOGq-2         ),N>& evk1,
				  const R_Q_square<2*(LOGQ-LOGq-2-  LOGDELTA),N>& evk2,
				  const R_Q_square<2*(LOGQ-LOGq-2-2*LOGDELTA),N>& evk3,
				  const R_Q_square<2*(LOGQ-LOGq-2-3*LOGDELTA),N>& evk4,
			      const R_Q_square<2*(LOGQ-LOGq-2-4*LOGDELTA),N>& evk5,
						R_Q_square<  (LOGQ-LOGq-2-5*LOGDELTA),N>& p    )
{	R_Q_square<LOGQ-LOGq-2,N> ct1;
	RS<LOGQ,LOGQ-LOGq-2,N>(ct,ct1); // Why rescale?
	ct1.print();
	const double u[32] = { 0.000000000000,0.140829338704,0.000000000000,-0.071956680803,
						   0.000000000000,0.314734688318,0.000000000000,-0.614678201819,
					      -0.000000000000,0.174769865184,-0.000000000000,0.037343013007,
					      -0.000000000000,-0.083190952031,-0.000000000000,0.128405049045,
						  -0.000000000000,0.247970475294,0.000000000000,-0.195489024287,
						   0.000000000000,0.784513028614,0.000000000000,-0.995925646437,
						   0.000000000000,0.336175543041,0.000000000000,-0.246866274711,
						 -0.000000000000,0.083032826731,-0.000000000000,-0.040430268293};
	eval_poly_deg31<LOGQ-LOGq-2,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,p);
	p.print();
	p *= 1ULL<<(LOGq+2);
}


template<int LOGQ, int N, int LOGq, int LOGDELTA>
void EvalMod_K2( const R_Q_square<  LOGQ,N>& ct,
				  const R_Q_square<2*(LOGQ-LOGq-1         ),N>& evk1,
				  const R_Q_square<2*(LOGQ-LOGq-1-  LOGDELTA),N>& evk2,
				  const R_Q_square<2*(LOGQ-LOGq-1-2*LOGDELTA),N>& evk3,
				  const R_Q_square<2*(LOGQ-LOGq-1-3*LOGDELTA),N>& evk4,
			      const R_Q_square<2*(LOGQ-LOGq-1-4*LOGDELTA),N>& evk5,
						R_Q_square<  (LOGQ-LOGq-1-5*LOGDELTA),N>& p    )
{	R_Q_square<LOGQ-LOGq-1,N> ct1;
	RS<LOGQ,LOGQ-LOGq-1,N>(ct,ct1);
	ct1.print();
	const double u[32] = {0.000000000000,0.132748607764,0.000000000000,-0.380278105603,0.000000000000,0.162806166803,-0.000000000000,0.107579389804,-0.000000000000,0.417564594173,-0.000000000000,-0.575577168110,-0.000000000000,0.279185145530,-0.000000000000,-0.149544709365,-0.000000000000,0.007450388247,0.000000000000,-0.002135232830,0.000000000000,0.000243436777,-0.000000000000,-0.000043966660,0.000000000000,0.000001649668,0.000000000000,-0.000000206220,-0.000000000000,0.000000011034,-0.000000000000,-0.000000001015 };
	eval_poly_deg31<LOGQ-LOGq-1,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,p);
	p.print();
	p *= 1ULL<<(LOGq);
}

template<int LOGQ, int N, int LOGq, int LOGDELTA>
void EvalMod_K3( const R_Q_square<  LOGQ,N>& ct,
				  const R_Q_square<2*(LOGQ-LOGq-  LOGDELTA),N>& evk1,
				  const R_Q_square<2*(LOGQ-LOGq-2*LOGDELTA),N>& evk2,
				  const R_Q_square<2*(LOGQ-LOGq-3*LOGDELTA),N>& evk3,
				  const R_Q_square<2*(LOGQ-LOGq-4*LOGDELTA),N>& evk4,
			      const R_Q_square<2*(LOGQ-LOGq-5*LOGDELTA),N>& evk5,
						R_Q_square<  (LOGQ-LOGq-6*LOGDELTA),N>& p    )
{	
	R_Q_square<  LOGQ,N> ct_copy=ct; ct_copy*=uint64_t((1ULL<<LOGDELTA)/3.);
	R_Q_square<LOGQ-LOGq-LOGDELTA,N> ct1;
	RS<LOGQ,LOGQ-LOGq-LOGDELTA,N>(ct_copy,ct1);
	ct1.print();
	const double u[32] = {-0.000000000000,0.137852010310,0.000000000000,-0.018904701355,0.000000000000,0.024320661898,-0.000000000000,-0.151590145511,-0.000000000000,0.367515735123,-0.000000000000,0.056870745879,-0.000000000000,0.541162271614,-0.000000000000,-1.048431384285,-0.000000000000,0.316696874787,0.000000000000,-0.282534155918,0.000000000000,0.106543202085,-0.000000000000,-0.053509305951,0.000000000000,0.005651171532,0.000000000000,-0.001841556587,-0.000000000000,0.000257314691,-0.000000000000,-0.000059365641};
	eval_poly_deg31<LOGQ-LOGq-LOGDELTA,N,LOGDELTA>(u,ct1,evk1,evk2,evk3,evk4,evk5,p);
	p.print();
	p *= 1ULL<<(LOGq);
}
*/
