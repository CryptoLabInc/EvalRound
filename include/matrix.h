#pragma once


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
			zero[s] = true;
			for (int i = 0; i < N; i++) vec[s][i] = 0;
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
					printf("%f", vec[s][i]);
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

