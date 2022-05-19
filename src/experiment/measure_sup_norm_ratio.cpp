#include <iostream>

#include "HEAAN/HEAAN.h"
#include "experiment/big.h"

void measure_message_sup_norm_ratio()
{
	for(int r = 0; r < 100; ++r) {
		Message<LOGN> z;
		set_random_message(z);
		double avg_norm = norm(z) * (2 / sqrt(N));
		double ratio = sup_norm(z) / avg_norm;
		std::cout << ratio << std::endl;
	}
}

void measure_plaintext_sup_norm_ratio()
{
	for(int r = 0; r < 100; ++r) {
		Message<LOGN> z;
		set_random_message(z);
		R_Q<LOGQ, N> pt;
		encode(z, Delta, pt);

		double norm = 0, sup_norm = 0;
		for(int i = 0; i < N; ++i) {
			Z_Q<LOGQ> val = pt[i];
			if(val.is_bigger_than_halfQ())
				val.negate();
			double val_abs_double = (double) val[0];
			norm += val_abs_double * val_abs_double;
			sup_norm = val_abs_double > sup_norm ? val_abs_double : sup_norm;
		}
		norm = sqrt(norm);

		double avg_norm = norm * (1 / sqrt(N));
		double ratio = sup_norm / avg_norm;
		std::cout << ratio << std::endl;
	}
}

void measure_norm_per_delta()
{
	for(int r = 0; r < 100; ++r) {
		Message<LOGN> z;
		set_random_message(z);
		R_Q<LOGQ, N> pt;
		encode(z, Delta, pt);

		double pt_norm = 0;
		double pt_sup_norm = 0;
		for(int i = 0; i < N; ++i) {
			Z_Q<LOGQ> val = pt[i];
			if(val.is_bigger_than_halfQ())
				val.negate();
			double val_abs_double = (double) val[0];

			pt_sup_norm = val_abs_double > pt_sup_norm ? val_abs_double : pt_sup_norm;
			pt_norm += val_abs_double * val_abs_double;
		}
		pt_norm = sqrt(pt_norm);

		//double ratio = (pt_sup_norm / Delta) / (sup_norm(z) / sqrt(N));
		//std::cout << ratio << std::endl;
		std::cout << (pt_norm / Delta) << std::endl;
	}
}

int main()
{
	measure_message_sup_norm_ratio();
	measure_plaintext_sup_norm_ratio();
	measure_norm_per_delta();
}