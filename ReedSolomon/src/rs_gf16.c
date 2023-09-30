
#include <stdint.h>
#include <stdio.h>
#include "rs_gf16.h"

#define PRIME4 0b10011	//the prime polynomial for GF(2^4), x^4 + x^1 + 1
#define SYMBOL_MASK4 0x0F0F0F0F0F0F0F0F	//masks off the low nibbles of each byte

const int8_t gf16_exp[] = {	// not a multiple of 2 so duplicate entries + offset needed for easy wraparound of negatives
	0x1,0x2,0x4,0x8,0x3,0x6,0xC,0xB,0x5,0xA,0x7,0xE,0xF,0xD,0x9,
	0x1,0x2,0x4,0x8,0x3,0x6,0xC,0xB,0x5,0xA,0x7,0xE,0xF,0xD,0x9
};

const int8_t gf16_log[] = {	// log_0 undefined so dummy 0xFF included to simplify indexing
	0xFF,0x0,0x1,0x4,0x2,0x8,0x5,0xA,0x3,0xE,0x9,0x7,0x6,0xD,0xB,0xC
};

int8_t gf16_div(int8_t a, int8_t b)
{
	if (b == 0)
		return -1;	// divide by 0 error, normal operation should never get here
	if (a == 0)
		return 0;
	
	return gf16_exp[gf16_log[a] - gf16_log[b] + 15];	// +15 offset to keep range positive
}

int8_t gf16_mul(int8_t a, int8_t b)
{
	if (a ==0 || b == 0)
		return 0;
	
	return gf16_exp[gf16_log[a] + gf16_log[b]];
}

int8_t gf16_pow(int8_t x, int8_t power)
{
	return gf16_exp[(gf16_log[x] * power) % 30];	// 30 is gf16_exp table size 
}

// slight optimization since most calls use x = 2 which evaluates to 1
int8_t gf16_2pow(int8_t power)
{	//TODO: check if it still needs the modulo for this special use case
	if (power >= 30)
		printf("required\n");
	return gf16_exp[power % 30];
}

int8_t gf16_inverse(int8_t x)
{
	return gf16_exp[15 - gf16_log[x]];
}

int8_t gf16_mul2_noLUT(int8_t x)
{
    x <<= 1;
    if (x > 15)  // if it's not within the bounds of the Galois Field, reduce by the primitive polynomial
        x ^= PRIME4;
    
    return x;
}

// helper for reducing after polynomial multiply
// TODO: check if inlining helps performance due to multiple issue/pipelining benefits
uint64_t gf16_poly_reduce(uint64_t p)
{
	uint64_t mask = p & 0xF0F0F0F0F0F0F0F0;
	mask |= (mask >> 3) ^ (mask >> 4);
	return p ^ mask;
}

/*
x86 has CLMUL for this type of operation and ARM NEON has VMUL in conjuction with P8 and P16 (aka polynomial)
types in newer architectures however it can be faked on other devices by creative packing of wider data types 
to acheieve a psuedo vectorized version that operates on 4 or 8 bytes in parallel depending on OS and hardware.

This works for GF(16) and GF(8) because polynomial multiply can be faked as a single long integer multiply so
long as there's sufficient padding between symbols allowing the overflow from one symbol to not clobber the next.
the required amount for polynomial scaling is 1 symbol width and for polynomial multiply is:
1 symbol width + log2(max result terms) - 1

Since indexing is only fast on byte bounds, we already have this padding for free and all that is needed is
a fast modular reduction which isn't too hard to encode as a series of shifts, masks, and xors but is unique
to any given primitive polynomial with diminishing returns for longer ones so usage of this on larger symbols
probably isn't practical.
*/
union Poly16 gf16_poly_mul(union Poly16 p, union Poly16 q)
{	
	p.uxl *= q.uxl;
	p.ull[0] = gf16_poly_reduce(p.ull[0]);
	p.ull[1] = gf16_poly_reduce(p.ull[1]);

	return p;
}

union Poly16 gf16_poly_scale(union Poly16 p, int8_t x)
{
	p.uxl *= x;
	
	p.ull[0] = gf16_poly_reduce(p.ull[0]);
	p.ull[1] = gf16_poly_reduce(p.ull[1]);

	return p;
}
