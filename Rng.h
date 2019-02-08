#if !defined __MY_RNG_H_
#define __MY_RNG_H_

// The RNG has been placed  outside of the main file. We do not
// 'cheat', and look at these variables, but I wanted to get the
// file completely out of the original source, to make it more
// apparent that we are FINDING these values, based upon several
// consecutive return values from the RNG

// even though this data is in code, we treat it like we do not know it.
// we want to use the above function to FIND these data.
uint64_t R_seed;
uint64_t RA = 6888241; // normally would be a const, BUT here we will modify them from within the check function.
uint64_t RK = 166671;
// other values to try for M.
//uint64_t RM = 4294967291;
//uint64_t RM = 0x100000001;		// care is needed if only 64 bit math used!!!
//uint64_t RM = 226888241;
//uint64_t RM = 65003;
//uint64_t RM = 0xFFFFFFFFFFFFFFFD;	// NOTE, will not work without 128 bit math
//uint64_t RM = 0xFFFFFFFFFFFFFFFF;
uint64_t RM = 0xFFFFFFFD;

// this is the 'real' PRNG function. We mess with the RA, RK, RM to use
// different values (in the test), set the seed, and simply call this
// function.  We then try to find RA, RK, RM solely based upon returned data.
uint64_t lcmprng() {
	// use 64 bit math for internal computations, to avoid overflow.
	uint64_t t = R_seed;
	R_seed = (t*RA + RK) % RM;
	return R_seed;
}

void lcmprng_srand(uint64_t v) {
	if (v == -1) {
		R_seed = time(0);
		return;
	}
	R_seed = v;
}


// here is the VC rand function.  We HAVE to figure this one out, on our quest to figure out
// the myspace random character data issues.

uint32_t vc_seed;

void msvc_srand(uint32_t x) {
	vc_seed = x;
}
int msvc_rand() {
	// here, m==2^32  a==214013 and k==2531011
	// but vc does a ((seed>>16)&0x7fff) to return only PART of the seed.
	return(((vc_seed = vc_seed * 214013L + 2531011L) >> 16) & 0x7fff);
}

#endif


