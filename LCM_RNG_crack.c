// This program will use a few consecutive returns from
// a linear congruent modular pRNG function, and it will
// find the 3 variables used:   The a, b, m.  The LCM
// is simply:   s1 = (s0*a+k)%M   s0 is the seed.  s1 is
// the next value returned (and seed for next random
// number).  This program will find the a, k and M
// it can do this with only 4 consecutive values.

// finds a, k and M.  See paper:  http://www.reteam.org/papers/e59.pdf

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>

#ifdef _MSC_VER
// my 64 bit GMP build against cygwin sux linking to 64 bit VC,
// but it 'can' be done, lol.
#pragma warning( disable : 4146)
#pragma warning( disable : 4244)
struct _reent * __getreent(void) {
	static char b[4096];
	return (struct _reent*)b;
}
#endif

#include <gmp.h>

int LCM_RNG_crack(uint64_t *v, int cnt, uint64_t *a, uint64_t *k, uint64_t *m);
static void mat_determ_3(uint64_t m[9], mpz_t ret);
static void gcd(mpz_t M, mpz_t A, mpz_t B);
static int mod_inv(mpz_t ret, mpz_t x, mpz_t M);

extern int dbg; // show more output.

// since GMP does not have mpz_set_ui64, we MAKE one NOTE, this works ONLY on 64 bit builds of GMP!
// it 'can' be done on 32 bit, but there, we do not assign to ._mp_d[0], bit instead typecast the pointer
// there to 64 bit, and assign to THAT pointer. Then on 32  bit, we have to set _mp_size to 1 or 2, depending
// upon if b is >= 2^32 or not.
// NOTE, if b==0, this function FAILS (I think).  I am not sure gmp works properly with size==1 and mp_d==0
// it 'may' handle that, but it may be undefined.
#define mpz_set_ui64(a,b) a[0]._mp_d[0] = b; a[0]._mp_size = 1


// This function finds a, k and M for lcm randome generator:  x1=(x0*a+k)%M; 
// it will first find M.  This is found using the determinate of a 3x3 matrix. That
// determinate should return M*x (x is some integer value).  We run 2 of these
// matrix computations, find the gcd of the 2 returns, and then factor out some
// tiny primes. Hopefully when done, we have the 'real' M value.  We then use that
// to find a and k, and check to make sure it produces proper data.  IF it does not
// produce proper data, we recursively call this function (if we can), throwing away
// the first v[0] item, so we get different 3x3 matrix and can hopefully find the
// correct M.
int LCM_RNG_crack(uint64_t *v, int cnt, uint64_t *_a, uint64_t *_k, uint64_t *_m) {
	mpz_t A, B, M, T, t1, t2, t3, t4, t5;
	int pr[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,71,73,79,83,89,97,101,103,107,109,113,127,131 };
	int i, j, reduced = 0, bad = 0;
	uint64_t M1[9] = { v[0], v[1], 1, v[1], v[2], 1, v[2], v[3], 1 };
	uint64_t M2[9] = { v[1], v[2], 1, v[2], v[3], 1, v[3], v[4], 1 };
	char Buf[128], Buf2[128];

	if (cnt < 5) {
		if (dbg)
			fprintf(stderr, "Error, solution NOT found!\n");
		return 0;
	}
	mpz_init_set_ui(A, 1); mpz_init_set_ui(B, 1); mpz_init_set_ui(M, 1); mpz_init_set_ui(T, 1);
	mpz_init_set_ui(t1, 1); mpz_init_set_ui(t2, 1); mpz_init_set_ui(t3, 1); mpz_init_set_ui(t4, 1);
	mpz_init_set_ui(t5, 1);

	// First off, find M.
	mat_determ_3(M1, A);
	mat_determ_3(M2, B);
	mpz_set(t1, A);
	mpz_set(t2, B);
	gcd(M, t1, t2);

	// note, M should be m*x where x is positive integer (could be 1)
	if (dbg>1) {
		mpz_get_str(Buf, 10, A); printf("dbg: matrix1 = A' = %s\n", Buf);
		mpz_get_str(Buf, 10, B); printf("dbg: matrix2 = B' = %s\n", Buf);
		mpz_get_str(Buf, 10, M); mpz_get_str(Buf2, 16, M); printf("dbg: M' = gcd(A',B')=%s  0x%s\n", Buf, Buf2);
	}
	// ok, try to reduce M, by factoring small primes
	// NOTE, if original m used in LCM was not prime (such as 0xFFFFFFFF),
	// then there can be small primes in it. We can NOT remove those!!
	// the check making sure M/pr[x] > than all incoming candidates is the
	// check which will hopefully catch these 'to be kept' small primes.
	// NOTE2, that check will fail sometimes, and small primes will be
	// removed, thus we may have to recursively run this function for win.
	for (i = (sizeof(pr) / sizeof(pr[0])) - 1; i >= 0; --i) {
		if (mpz_tdiv_ui(M, pr[i]) == 0) {
			bad = 0;
			mpz_divexact_ui(T, M, pr[i]);
			for (j = 0; j < cnt && !bad; ++j) {
				if (mpz_cmp_ui(T, v[j]) < 0)
					bad = 1;
			}
			if (!bad) {
				mpz_divexact_ui(M, M, pr[i]);
				if (dbg>1)
					printf("dbg: removed factor %d\n", pr[i]);
				++i;
				reduced = 1;
			}
		}
	}

	if (M[0]._mp_size > 1) {
		mpz_clear(A); mpz_clear(B); mpz_clear(M); mpz_clear(T);
		mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4); mpz_clear(t5);
		if (dbg>1)
			printf("dbg: The M value found is > 64 bits.  Assuming failure, and trying again\n");
		if (cnt <= 5)
			if (cnt < 5) {
				if (dbg)
					fprintf(stderr, "Error, solution NOT found with data provided! (failed size of m > 64 bits)\n");
				return 0;
			}
		return LCM_RNG_crack(&v[1], cnt - 1, _a, _k, _m);
	}

	*_m = M[0]._mp_d[0];
	if (dbg)
		printf("m = %"PRIu64"  0x%"PRIx64"\n", *_m, *_m);

	if (dbg>1) {
		printf("dbg: Finding a and k, then validating that the lcm-prng works properly using for them. If\n");
		printf("dbg: they do not work, we recursively call lcmprng_crack until out of values, solution found\n");
		printf("dbg: prn[0] = %"PRIu64"\n", v[0]);
		printf("dbg: prn[1] = %"PRIu64"\n", v[1]);
		printf("dbg: prn[2] = %"PRIu64"\n", v[2]);
		printf("dbg: to find a and k, we solve these 2 modular formula\n");
		printf("dbg: v[0]a+k == v[1] mod(m)   v[1].a+k == v[2] mod(m)   solve for a and k\n");
		printf("dbg: %"PRIu64" == %"PRIu64".a+k mod(%s)\n", v[1], v[0], Buf);
		printf("dbg: %"PRIu64" == %"PRIu64".a+k mod(%s)\n", v[2], v[1], Buf);
	}
	// now find a and k
	// EX1 = v[1] == v[0]a+k (mod m)   EX2 = v[2] == v[1].a+k (mod m)   solve for a and k
	// subtract EX1 - EX2 (mod m)
	mpz_set_ui64(t1, v[1]);
	mpz_set_ui64(t2, v[2]);
	mpz_sub(t1, t1, t2);
	if (mpz_cmp_ui(t1, 0) < 0) mpz_neg(t1, t1);
	mpz_set_ui64(t2, v[0]);
	mpz_set_ui64(t3, v[1]);
	mpz_sub(t2, t2, t3);
	if (mpz_cmp_ui(t2, 0) < 0) mpz_neg(t2, t2);
	// now t1 == t2.a (mod m)  (note, k is now gone!)
	// we now find  t2^-1 mod m, to get 1.a (i.e. a = result) again, all (mod m)
	// this is multiplicative inverse, we find it we e-euclidean algorithm.
	if (!mod_inv(t3, t2, M)) {	// now, t3 == t2^-1 (mod m)
		mpz_clear(A); mpz_clear(B); mpz_clear(M); mpz_clear(T);
		mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4); mpz_clear(t5);
		if (cnt <= 5) {
			if (dbg)
				fprintf(stderr, "Error, solution NOT found with data provided due to inability to find mod_inverse!\n");
			return 0;
		}
		return LCM_RNG_crack(&v[1], cnt - 1, _a, _k, _m);
	}
	mpz_mul(t4, t1, t3);	// compute t1 * t2^-1 which is the VALUE of a
	mpz_mod(t5, t4, M);
	*_a = t5[0]._mp_d[0];
	if (dbg)
		printf("a = %"PRIu64"  0x%"PRIx64"\n", *_a, *_a);

	// now plug a into the 2nd formula, EX2.  i..e k = v[2]-v[1]*a  (mod m)
	mpz_set_ui64(t1, v[2]);
	mpz_set_ui64(t2, v[1]);
	mpz_set_ui64(t3, *_a);
	mpz_mul(t2, t2, t3);
	mpz_sub(t1, t1, t2);
	mpz_mod(t1, t1, M);
	if (mpz_cmp_ui(t1, 0) < 0) mpz_neg(t1, t1);
	*_k = t1[0]._mp_d[0];
	if (dbg)
		printf("k = %"PRIu64"  0x%"PRIx64"\n", *_k, *_k);

	// now make sure that a,k,m 'work', by using it against the passed in data.
	{
		uint64_t seed = v[0];
		int i;
		// NOTE, we use GMP here, since m can be up to 2^64-1, and 
		// we must avoid overflow problems.
		mpz_set_ui64(t2, seed);
		mpz_set_ui64(t3, *_a);
		mpz_set_ui64(t4, *_k);
		mpz_set_ui64(t5, *_m);
		for (i = 1; i < cnt; ++i) {
			//seed = (seed * *_a + *_k) % *_m;
			mpz_mul(t1, t2, t3);
			mpz_add(t1, t1, t4);
			mpz_mod(t2, t1, t5);
			seed = t2[0]._mp_d[0];
			if (seed != v[i]) {
				mpz_clear(A); mpz_clear(B); mpz_clear(M); mpz_clear(T);
				mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4); mpz_clear(t5);
				if (dbg>1)
					printf("FAILED, at index %i\n", i);
				if (cnt <= 5)
					if (cnt <= 5) {
						if (dbg)
							fprintf(stderr, "Error, solution NOT found with data provided (failed in playback)!\n");
						return 0;
					}
				return LCM_RNG_crack(&v[1], cnt - 1, _a, _k, _m);
			}
		}
	}
	printf("\n");
	printf("******LCM-PRNG cracked******\n");
	printf(" Function: nx[n+1] = (x[n]*a + k) %% m;  where\n");
	printf("   m=%"PRId64"\n   a=%"PRId64"\n   k=%"PRId64"\n", *_m, *_a, *_k);

	mpz_clear(A); mpz_clear(B); mpz_clear(M); mpz_clear(T);
	mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4); mpz_clear(t5);

	return 1;
}

static void gcdExtended(mpz_t g, const mpz_t a, const mpz_t b, mpz_t x, mpz_t y)
{
	mpz_t x1, y1, g1, b1, t1; // To store results of recursive call 
	// Base Case 
	if (mpz_cmp_ui(a, 0) == 0) {
		mpz_set_ui(x, 0);
		mpz_set_ui(y, 1);
		mpz_set(g, b);
		return;
	}
	mpz_init_set_ui(x1, 1); mpz_init_set_ui(y1, 1); mpz_init_set_ui(g1, 1); mpz_init_set_ui(b1, 1); mpz_init_set_ui(t1, 1);

	mpz_mod(b1, b, a);
	gcdExtended(g1, b1, a, x1, y1);

	// Update x and y using results of recursive 
	// call 
	//*x = y1 - (b / a) * x1;
	mpz_div(t1, b, a);
	mpz_mul(t1, t1, x1);
	mpz_sub(x, y1, t1);
	//*y = x1;
	mpz_set(y, x1);

	//return gcd;
	mpz_set(g, g1);
	mpz_clear(x1); mpz_clear(y1); mpz_clear(g1); mpz_clear(b1); mpz_clear(t1);
}

static int mod_inv(mpz_t ret, mpz_t a, mpz_t M) {
	//mpz_set_ui64(ret, 3374592081);
	//   3374592081 known to work for lcmprng_srand(888); 
	//   ONLY used until full mod_inv function was completed
	mpz_t x, y, g;

	mpz_init_set_ui(x, 1); mpz_init_set_ui(y, 1); mpz_init_set_ui(g, 1);
	gcdExtended(g, a, M, x, y);
	if (mpz_cmp_ui(g, 1) != 0) {
		if (dbg>1)
			fprintf(stderr, "Error, mod_inv FAILS!\n");
		mpz_clear(x); mpz_clear(y); mpz_clear(g);
		return 0;
	}
	// make sure ret (the inverse), is >= 0 and < M, i.e. x=(x%m+m)%m
	mpz_mod(x, x, M);
	mpz_add(x, x, M);
	mpz_mod(ret, x, M);
	mpz_clear(x); mpz_clear(y); mpz_clear(g);
	return 1;
}
static void gcd(mpz_t M, mpz_t n1, mpz_t n2) {
	if (mpz_cmp_ui(n1, 0) < 0)
		mpz_neg(n1, n1);
	if (mpz_cmp_ui(n2, 0) < 0)
		mpz_neg(n2, n2);
	while (mpz_cmp(n1, n2))
	{
		if (mpz_cmp(n1, n2) > 0)
			mpz_sub(n1, n1, n2);
		else
			mpz_sub(n2, n2, n1);
	}
	mpz_set(M, n1);
}
/*
*      +-     -+
*      | a b c |
*  A = | d e f |
*      | g h i |
*      +-     -+
*
*  A = a(ei-fh) - b(di-fg) + c(dh-eg)
*/
static void mat_determ_3(uint64_t m[9], mpz_t A) {
	uint64_t a, b, c, d, e, f, g, h, i;
	mpz_t t1, t2, t3, t4;
	a = m[0]; b = m[1]; c = m[2];
	d = m[3]; e = m[4]; f = m[5];
	g = m[6]; h = m[7]; i = m[8];

	mpz_init_set_ui(t1, 1); mpz_init_set_ui(t2, 1); mpz_init_set_ui(t3, 1); mpz_init_set_ui(t4, 1);

	// set A to : a(ei-fh)
	mpz_set_ui64(t1, a);	// t1 = a
	mpz_set_ui64(t2, e);
	mpz_set_ui64(t3, i);
	mpz_mul(t2, t2, t3);	// t2 = ei
	mpz_set_ui64(t3, f);
	mpz_set_ui64(t4, h);
	mpz_mul(t3, t3, t4);	// t3 = fh
	mpz_sub(t2, t2, t3);    // t2 = ei-fh
	mpz_mul(A, t1, t2);	// A = a(ei-fh)

				// A -= b(di-fg)
	mpz_set_ui64(t1, b);	// t1 = b
	mpz_set_ui64(t2, d);
	mpz_set_ui64(t3, i);
	mpz_mul(t2, t2, t3);	// t2 = di
	mpz_set_ui64(t3, f);
	mpz_set_ui64(t4, g);
	mpz_mul(t3, t3, t4);	// t3 = fg
	mpz_sub(t2, t2, t3);    // t2 = di-fg
	mpz_mul(t1, t1, t2);	// t1 = b(di-fg)
	mpz_sub(A, A, t1);	// A = a(ei-fh) - b(di-fg)

				// A += c(dh-eg)
	mpz_set_ui64(t1, c);	// t1 = c
	mpz_set_ui64(t2, d);
	mpz_set_ui64(t3, h);
	mpz_mul(t2, t2, t3);	// t2 = dh
	mpz_set_ui64(t3, e);
	mpz_set_ui64(t4, g);
	mpz_mul(t3, t3, t4);	// t3 = eg
	mpz_sub(t2, t2, t3);    // t2 = dh-eg
	mpz_mul(t1, t1, t2);	// t1 = c(dh-eg)
	mpz_add(A, A, t1);	// A = a(ei-fh) - b(di-fg) + c(dh-eg)

	mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4);

	// A == a(ei-fh) - b(di-fg) + c(dh-eg)
	// A is the return value.
}