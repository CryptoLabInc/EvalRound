#include "HEAAN/HEAAN_matrix.h"
#include "HEAAN/DFT.h"

int main() {
	SparseDiagonal<1<<9,3> Vr[9];
	SparseDiagonal<1<<9,3> Vi[9];
	splitU0NR<10>(Vr,Vi);
	SparseDiagonal<1<<9,9> A0r,A0i;
	SparseDiagonal<1<<9,9> A1r,A1i;
	SparseDiagonal<1<<9,9> A2r,A2i;
	SparseDiagonal<1<<9,9> tempr,tempi;
	SparseDiagonal<1<<9,27> A3r,A3i;
	MatMul<1<<9,3,3>(Vr[0],Vi[0],Vr[1],Vi[1],A0r,A0i);
	MatMul<1<<9,3,3>(Vr[2],Vi[2],Vr[3],Vi[3],A1r,A1i);
	MatMul<1<<9,3,3>(Vr[4],Vi[4],Vr[5],Vi[5],A2r,A2i);
	MatMul<1<<9,3,3>(Vr[6],Vi[6],Vr[7],Vi[7],tempr,tempi);
	MatMul<1<<9,9,3>(tempr,tempi,Vr[8],Vi[8],A3r,A3i);
	A0r.print("A0r");A0i.print("A0i");
	A1r.print("A1r");A1i.print("A1i");
	A2r.print("A2r");A2i.print("A2i");
	A3r.print("A3r");A3i.print("A3i");
}