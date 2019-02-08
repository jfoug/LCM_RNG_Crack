// main() for LCM_RNG random number 'cracker'
//
// simple LCM is:    x[n+1] = (x[n]*a+k)%M;
//
// This program finds a, k and M, based only on
// a set of consecutive x[] values
//
// See paper:  http://www.reteam.org/papers/e59.pdf
////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include "Rng.h"

int LCM_RNG_crack(uint64_t *v, int cnt, uint64_t *a, uint64_t *k, uint64_t *m);

int dbg = 0; // show more output.

int main(int argc, char **argv) {
	uint64_t v[24], a, k, m;
	int i, j;

	//lcmprng_srand(-1); // seed to an unknown value;
	lcmprng_srand(889);
	// ok, get some data from the RNG
	for (i = 0; i < 24; ++i)
		v[i] = lcmprng();
	// now try to crack the RNG
	if (!LCM_RNG_crack(v, 24, &a, &k, &m))
		printf("We could not crack the LCM!\n");

	// ok, try to now crack VC's random generator
	srand(888);
	msvc_srand(888);
	for (i = 0; i < 24; ++i) {
		uint32_t x = msvc_rand();
		v[i] = rand();
		if (x != v[i])
			printf("");
	}

	printf("\n\nNow trying to crack VC's rand function\n");
	if (!LCM_RNG_crack(v, 24, &a, &k, &m))
		printf("We could not crack the LCM!\n");

	printf("\nDone with test\n\n");
}

