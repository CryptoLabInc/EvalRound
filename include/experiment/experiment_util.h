#include "impl/message.h"
#include "HEAAN/DFT.h"

template <int LOGN>
void set_test_message(Message<LOGN> &z);

template <int LOGN>
void set_random_message(Message<LOGN> &z);

template <int LOGN>
void set_test_rounded_message(Message<LOGN> &z);

template <int LOGN, int K>
void set_test_matrix(Message<LOGN> &A[K]);

template <int LOGN, int K>
void set_random_matrix(Message<LOGN> &A[K]);

template <int LOGN>
void set_test_U0_matrix(SparseDiagonal<(1<<(LOGN-1)),3> U0r[LOGN-1],
	            SparseDiagonal<(1<<(LOGN-1)),3> U0i[LOGN-1]);