#include "gf8.h"
#include <stdio.h>  //TODO: remove this and the associated printf
#define PRIME_GF8 0b1011	//the prime polynomial for GF(8), x^3 + x^1 + 1

//simplified galois field multiply by 2 used for generating the Look Up Tables
int8_t gf8_mul2_noLUT(int8_t x)
{
    x <<= 1;
    if (x > 7)  // if it's not within the bounds of the Galois Field, reduce by the primitive polynomial
        x ^= PRIME_GF8;
    
    return x;
}

int8_t gf8_div(int8_t a, int8_t b)
{
	if (b == 0)
		return -1;	// divide by 0 error, normal operation should never get here
	if (a == 0)
		return 0;
	
	return gf8_exp[gf8_log[a] - gf8_log[b] + 7];	// +7 offset to keep range positive
}

int8_t gf8_mul(int8_t a, int8_t b)
{
	if (a == 0 || b == 0)
		return 0;
	
	return gf8_exp[gf8_log[a] + gf8_log[b]];
}

int8_t gf8_pow(int8_t x, int8_t power)
{
	return gf8_exp[(gf8_log[x] * power) % 14];	// 14 is gf8_exp table size 
}

// slight optimization since most calls use x = 2 which evaluates to 1
int8_t gf8_2pow(int8_t power)
{	//TODO: check if it still needs the modulo for this special use case
	if (power >= 14)
		printf("required\n");
	return gf8_exp[power % 14];
}

int8_t gf8_inverse(int8_t x)
{
	return gf8_exp[15 - gf8_log[x]];
}

//prior to reduction, term can extend up to 2 bits above symbol due to shifting
//this function is customized to GF(8) with prime polynomial 1011
uint32_t gf8_poly_reduce(uint32_t p, uint32_t of)
{
	return p ^ (of >> 2) ^ (of >> 1);
}

// optimized for fewer memory accesses
//TODO: check if multiplies are faster, currently assuming that single shifts and conditional assignment are better
uint32_t gf8_poly_scale(uint32_t p, int8_t x)
{
	uint32_t r0, r1, r2, of;
	r0 = (x & 1) ? p : 0;
	p <<= 1;
	r1 = (x & 2) ? p : 0;
	p <<= 1;
	r2 = (x & 4) ? p : 0;

	of = (r1 & 011111111110) ^ (r2 & 033333333330);
	r0 ^= (r1 & 006666666666) ^ (r2 & 004444444444);

	return gf8_poly_reduce(r0, of);
}

//TODO: the resulting polynomial is longer than the inputs, might have to rethink the types if this exceeds 7 terms
// also might remove the size args since the loop continue should handle it fine and it's only a max array length of 7 for both

//Assumes that result can never be longer than 7 terms, and the shorter polynomial is in q
// currently assuming the second multiplier is no more than 3 terms, this is probably not enough
// may have to consider rearranging as interleaved high and low degree symbols if it results in faster ops
uint32_t gf8_poly_mul(uint32_t p, uint32_t q)
{
	uint32_t r0, r1, r2, of;
	r0 = (q & 01) ? p : 0;
	r1 = (q & 02) * p;
	r2 = (q & 04) * p;

	r0 ^= (q & 010) * p;
	r1 ^= (q & 020) * p;
	r2 ^= (q & 040) * p;

	r0 ^= (q & 0100) * p;
	r1 ^= (q & 0200) * p;
	r2 ^= (q & 0400) * p;

	of = (r1 & 011111111110) ^ (r2 & 033333333330);
	r0 ^= (r1 & 006666666666) ^ (r2 & 004444444444);

	return gf8_poly_reduce(r0, of);
}

//p is dividend, q is divisor, p_sz and q_sz are size in BITS not symbols
//returns remainder of the division since the quotient is never used
uint32_t gf8_poly_mod(uint32_t p, int8_t p_sz, uint32_t q, int8_t q_sz)
{
	//if p_sz and q_sz is known at compile time, this can be rewritten to be unrollable
	p_sz -= GF8_IDX_INC;
	q_sz -= GF8_IDX_INC;
	//uncomment the following line to return the quotient and remainder in a single return value with the start of the quotient at b_arr[q_sz - 2]
	//q &= ~(-1 << q_sz); //clears the highest order term which should be a 1
	p <<= q_sz;
	q <<= p_sz;
	for(int8_t i = p_sz + q_sz; i >= q_sz; i -= GF8_IDX_INC)
	{
		p ^= gf8_poly_scale(q, (p >> i) & 7);
		q >>= 3;
	}
	
	return p;
}

//optimized version of div for binomial divisor/single eval point
//TODO: check if this is actually more efficient at this size
int8_t gf8_poly_eval(uint32_t p, int8_t p_sz, int8_t x)
{
	p_sz -= GF8_IDX_INC;
    int8_t y = p >> p_sz;
	int8_t logx = gf8_log[x];
    for(p_sz -= GF8_IDX_INC; p_sz >= 0; p_sz -= GF8_IDX_INC)
        y = gf8_exp[gf8_log[y] + logx] ^ ((p >> p_sz) & 7);

    return y;
}