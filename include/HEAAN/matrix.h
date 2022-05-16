#pragma once

#include <cstdio>
#include <cassert>
#include <iostream>

template<int N>
void matrix_product(const double Ar[N][N],
	                const double Ai[N][N],
	                const double Br[N][N],
	                const double Bi[N][N], double Cr[N][N],
	                                       double Ci[N][N]) {
	for (int i = 0; i<N; i++)
	for (int j = 0; j<N; j++) {
		Cr[i][j] = 0;
		Ci[i][j] = 0;
		for (int k = 0; k<N; k++) {
			Cr[i][j] += Ar[i][k] * Br[k][j] - Ai[i][k] * Bi[k][j];
			Ci[i][j] += Ar[i][k] * Bi[k][j] + Ai[i][k] * Br[k][j];
		}
	}
}

template<int N>
void print(const double Ar[N][N],
	       const double Ai[N][N]) {
	printf("matrix(");
	for (int i = 0; i<N; i++) {
		printf("[");
		for (int j = 0; j<N; j++)
			if (j<N - 1) printf("%2.12f+%%i*%2.12f,", Ar[i][j], Ai[i][j]);
			else         printf("%2.12f+%%i*%2.12f]", Ar[i][j], Ai[i][j]);
			if (i<N - 1) printf(", \n");
			else         printf(");\n");
	}
}


//------------------------------------------------------------------------------------
// sparse diagonal matrix of size NxN
// with S number of diagonal vectors
//------------------------------------------------------------------------------------
template<int N, int S>
struct SparseDiagonal {
	double vec[S][N];
	int    off[S];
	bool  zero[S];

	// copy operator
	void operator=(const SparseDiagonal<N,S> &other) {
		for(int s = 0; s < S; ++s) {
			for(int i = 0; i < N; ++i) {
				vec[s][i] = other.vec[s][i];
			}
			off[s] = other.off[s];
			zero[s] = other.zero[s];
		}
	}
	//
	double get_element(int i, int j) const {
		for (int s = 0; s < S; s++)
			if ((j - i - off[s]) % N == 0)
				return vec[s][i];
		return 0;
	}
	//
	void print(const char* name) const {
		printf("%s:zeromatrix(%d,%d)$\n", name, N, N);
		for (int s = 0; s < S; s++)
		for (int i = 0; i < N; i++) {
				int j = (i + off[s]) % N;
				if (j < 0) j += N;
				if (vec[s][i] != 0)
					printf("%s[%d][%d]:%f$\n", name, i + 1, j + 1, vec[s][i]);
		}
	}
	//
	void setzero() {
		for (int s = 0; s < S; s++) {
			for (int i = 0; i < N; i++) vec[s][i] = 0;
			zero[s] = true;
			off[s] = 0; // for consistency
		}
	}
	//
	void set_element(int i, int j, double v) {
		for (int s = 0; s < S; s++)
			if ((j - i - off[s]) % N == 0) {
				vec[s][i] = v;
				zero[s] = false;
				return;
			}
		assert(false);
	}
	//
	void operator+=(const SparseDiagonal<N, S>& B) {
		for (int s = 0; s < S; s++) {
			assert(off[s] == B.off[s]);
			for (int i = 0; i < N; i++)
				vec[s][i] += B.vec[s][i];
		}
	}
	//
	void operator-=(const SparseDiagonal<N, S>& B) {
		for (int s = 0; s < S; s++) {
			assert(off[s] == B.off[s]);
			for (int i = 0; i < N; i++)
				vec[s][i] -= B.vec[s][i];
		}
	}
	//
	void operator*=(double a) {
		for (int s = 0; s < S; s++)
			for (int i = 0; i < N; i++)
				vec[s][i] *= a;
	}

	void transpose() {
		for(int s=0; s<S; s++)
			if (off[s] > 0) {
				double temp[N];
				for (int i = 0; i < N; i++) {
					if (i < N - off[s]) temp[off[s] + i] = vec[s][i];
					else				temp[i - (N - off[s])] = vec[s][i];
				}
				for (int i = 0; i < N; i++) vec[s][i] = temp[i];
				off[s] = N - off[s];
					
			}
	}

	void negate() {
		for (int s = 0; s < S; s++)
			for (int i = 0; i < N; i++)
				vec[s][i] = -vec[s][i];
	}
};

//------------------------------------------------------------------------------------------
// Matrix vector product of complex SparseDiagonal matrix and a complex vector
//------------------------------------------------------------------------------------------
template<int LOGN>
void matrix_vector_product(
    const double zr[1 << (LOGN - 1)], const double zi[1 << (LOGN - 1)],
    SparseDiagonal<(1<<(LOGN-1)),3> Ar,
	SparseDiagonal<(1<<(LOGN-1)),3> Ai,
    double Azr[1 << (LOGN - 1)], double Azi[1 << (LOGN - 1)]) {
	const int N = 1 << LOGN;
    for(int i = 0; i < N/2; ++i) {
        double sumr = 0, sumi = 0;
        for(int k = 0; k < 3; ++k) {
            int jr = (i+Ar.off[k]) % (N/2);
            int ji = (i+Ai.off[k]) % (N/2);
            sumr += Ar.vec[k][i] * zr[jr] - Ai.vec[k][i] * zi[ji];
            sumi += Ar.vec[k][i] * zi[ji] + Ai.vec[k][i] * zr[jr];
        }
        Azr[i] = sumr;
        Azi[i] = sumi;
    }
}


//------------------------------------------------------------------------------------------
// Matrix multiplication of two SparseDiagonal matrices
//------------------------------------------------------------------------------------------
template< int N, int S1, int S2 >
void MatMul(const SparseDiagonal<N, S1   >& A,
	        const SparseDiagonal<N,    S2>& B,
	              SparseDiagonal<N, S1*S2>& C) {
	int count = 0; C.setzero();
	for (int s1 = 0; s1 < S1; s1++)
	for (int s2 = 0; s2 < S2; s2++) {
		int off1 = A.off[s1];
		int off2 = B.off[s2];
		int off3 = (off1 + off2) % N;
		int s3 = -1;
		for (int k = 0; (k < count) && (s3 == -1); k++)
			if (C.off[k] == off3)
				s3 = k;
		if (s3 == -1) {
			s3 = count;
			C.zero[s3] = false;
			C.off [s3] =  off3;
			count++;
		}
		for (int i = 0; i < N; i++)
			C.vec[s3][i] += A.vec[s1][i] * B.vec[s2][(i + off1) % N];
	}
}

template< int N, int S1, int S2 >
void MatMul(const SparseDiagonal<N, S1   >&Ar,const SparseDiagonal<N, S1   >&Ai,
			const SparseDiagonal<N,    S2>&Br,const SparseDiagonal<N,    S2>&Bi,
				  SparseDiagonal<N, S1*S2>&Cr,      SparseDiagonal<N, S1*S2>&Ci){
	int count = 0;
	Cr.setzero(); Ci.setzero();
	for (int s1 = 0; s1 < S1; s1++)
	for (int s2 = 0; s2 < S2; s2++) {
		if(Ar.zero[s1] && Ai.zero[s1] || Br.zero[s2] && Bi.zero[s2])
			continue;
		assert(Ar.off[s1] == Ai.off[s1]);
		assert(Br.off[s2] == Bi.off[s2]);
		int off1 = Ar.off[s1];
		int off2 = Br.off[s2];
		int off3 = (off1 + off2) % N;
		int s3 = -1;
		for (int k = 0; (k < count) && (s3 == -1); k++) {
			assert(Cr.off[k] == Ci.off[k]);
			if (Cr.off[k] == off3)
				s3 = k;
		}
		if (s3 == -1) {
			s3 = count;
			Cr.zero[s3] = false; Ci.zero[s3] = false;
			Cr.off [s3] =  off3; Ci.off[s3] = off3;
			count++;
		}
		for (int i = 0; i < N; i++) {
			Cr.vec[s3][i] += Ar.vec[s1][i] * Br.vec[s2][(i + off1) % N];
			Cr.vec[s3][i] -= Ai.vec[s1][i] * Bi.vec[s2][(i + off1) % N];
			Ci.vec[s3][i] += Ar.vec[s1][i] * Bi.vec[s2][(i + off1) % N];
			Ci.vec[s3][i] += Ai.vec[s1][i] * Br.vec[s2][(i + off1) % N];
		}
	}
}

//------------------------------------------------------------------------------------------
// Remove unnecessary zero vectors
//------------------------------------------------------------------------------------------
template< int N, int S1, int S2 >
void RemoveClutter(const SparseDiagonal<N, S1>& A,
	                     SparseDiagonal<N, S2>& B) {
	assert(S2 < S1);
	B.setzero();
	for (int s = 0; s < S2; s++) {
		B.off[s] = A.off[s];
		B.zero[s] = false;
		for (int i = 0; i < N; i++)
			B.vec[s][i] = A.vec[s][i];
	}
	for (int s = S2; s < S1; s++)
		assert(A.zero[s] == true);
}

template <int N, int S>
void print_max(const SparseDiagonal<N, S>& A) {
	double max = 0;
	for (int s = 0; s < S; s++) {
		for (int i = 0; i < N; i++) {
			double a = std::abs(A.vec[s][i]);
			max = max > a ? max : a;
		}
	}
	std::cout << "Max " << max << std::endl;

}