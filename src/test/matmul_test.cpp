#include "experiment/rns_debug.h"

#include "HEAAN/HEAAN.h"
#include "HEAAN/HEAAN_matrix.h"

bool isEqual(double a, double b) {
	if(a == 0)
		return a == b;
	double diff_rate = (b - a) / a;
	return -0.001 < diff_rate && diff_rate < 0.001;
}

void element_test() {
	SparseDiagonal<1<<(LOGN-1), 3> U0r[LOGN-1];
	SparseDiagonal<1<<(LOGN-1), 3> U0i[LOGN-1];
	set_test_U0_matrix<LOGN>(U0r, U0i);

	SparseDiagonal<1<<(LOGN-1),9> Ar,Ai;
	MatMul<1<<(LOGN-1),3,3>(U0r[0],U0i[0],U0r[1],U0i[1],Ar,Ai);

	for(int i = 0; i < N/2; ++i) {
		for(int j = 0; j < N/2; ++j) {
			double sumr = 0, sumi = 0;
			for(int k = 0; k < N/2; ++k) {
				sumr += U0r[0].get_element(i, k) * U0r[1].get_element(k, j);
				sumr -= U0i[0].get_element(i, k) * U0i[1].get_element(k, j);
				sumi += U0r[0].get_element(i, k) * U0i[1].get_element(k, j);
				sumi += U0i[0].get_element(i, k) * U0r[1].get_element(k, j);
			}

			if(!isEqual(Ar.get_element(i, j), sumr) || !isEqual(Ai.get_element(i, j), sumi)) {
				std::cout << i << " " << j << std::endl;
				std::cout << Ar.get_element(i, j) << " " << sumr << std::endl;
				std::cout << Ai.get_element(i, j) << " " << sumr << std::endl;
			}
		}
	}
}

void vector_test() {
	SparseDiagonal<1<<(LOGN-1), 3> U0r[LOGN-1];
	SparseDiagonal<1<<(LOGN-1), 3> U0i[LOGN-1];
	set_test_U0_matrix<LOGN>(U0r, U0i);

	SparseDiagonal<1<<(LOGN-1),9> Ar,Ai;
	MatMul<1<<(LOGN-1),3,3>(U0r[1],U0i[1],U0r[0],U0i[0],Ar,Ai);

	Message<LOGN> z, U0z[2], Az;
	set_random_message(z);
	matrix_vector_product(z, U0r[0], U0i[0], U0z[0]);
	matrix_vector_product(U0z[0], U0r[1], U0i[1], U0z[1]);
	matrix_vector_product(z, Ar, Ai, Az);
	print("U0z", U0z[1]);
	print("Az", Az);
}

void linear_transform_test() {
	SparseDiagonal<1<<(LOGN-1), 3> U0r[LOGN-1];
	SparseDiagonal<1<<(LOGN-1), 3> U0i[LOGN-1];
	set_test_U0_matrix<LOGN>(U0r, U0i);

	int s[N];
    HEAAN<LOGQ,N>::keygen(H,s);
   
    Message<LOGN> z, Az;
    set_random_message(z);

    R_Q<LOGQ, N> pt;
    R_Q_square<LOGQ,N> ct, ct_ser, ct_gr2, ct_gr4;
    encode(z,Delta,pt);
    HEAAN<LOGQ,N>::enc(pt,s,ct);
	serial_linear_transform<LOGQ, LOGN, LOGDELTA_TILDE, 3, 8>(U0r, U0i, ct, s, ct_ser);
	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA_TILDE, 3, 8, 2>(U0r, U0i, ct, s, ct_gr2);
	grouped_serial_linear_transform<LOGQ, LOGN, LOGDELTA_TILDE, 3, 8, 4>(U0r, U0i, ct, s, ct_gr4);

	R_Q<LOGQ, N> pt_ser, pt_gr2, pt_gr4;
	Message<LOGN> z_ser, z_gr2, z_gr4;
    HEAAN<LOGQ,N>::dec(ct_ser, s, pt_ser);
	decode_log(pt_ser, LOGDELTA+8*LOGDELTA_TILDE, z_ser);
	HEAAN<LOGQ,N>::dec(ct_gr2, s, pt_gr2);
    decode_log(pt_gr2,LOGDELTA+4*LOGDELTA_TILDE,z_gr2);
	HEAAN<LOGQ,N>::dec(ct_gr4, s, pt_gr4);
    decode_log(pt_gr4,LOGDELTA+2*LOGDELTA_TILDE,z_gr4);
	print("serial linear transformed", z_ser);
	print("grouped by 2 serial linear transformed", z_gr2);
	print("grouped by 4 serial linear transformed", z_gr4);
}

int main() {
	linear_transform_test();	
}