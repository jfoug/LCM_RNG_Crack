default: linmod_rev.c LCM_RNG_crack.c Rng.h
	gcc -o linmod_rev linmod_rev.c LCM_RNG_crack.c -lgmp
