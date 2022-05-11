#include <iostream>

#include "include/HEAAN/advanced/mod.h"

int main()
{
	uint64_t a = 12473026286576736102ULL;
	uint64_t b = 4953360856192166699ULL;
	uint64_t lo_expected = 9994061768550322210ULL;
	uint64_t hi_expected = 3349284834185960393ULL;
	uint64_t lo = 0, hi = 0;
	mul(a, b, lo, hi);
	std::cout << lo - lo_expected << std::endl;
	std::cout << hi - hi_expected << std::endl;
}
