#include <iostream>

#include "HEAAN/advanced/mod.h"

int main()
{
	uint64_t al = 2487277757547979106ULL;
	uint64_t ah = 16209862899231289696ULL;
	uint64_t q = 8471933682407887875ULL;
	uint64_t expected = 3357710663169374592ULL;
	std::cout << expected - mod(al, ah, q) << std::endl;
}
