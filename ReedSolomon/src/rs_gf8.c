//BCH view, systematic encoding Reed Solomon using 3 bit symbols
#include <stdint.h>
#include <stdio.h>
#include "rs_gf8.h"

#define PRIME_GF8 			0b1011	//the prime polynomial for GF(8), x^3 + x^1 + 1
#define OCTAL_MASK_A 		00707070707
#define OCTAL_MASK_B 		07070707070
#define GF8_IDX_INC			3	//how many bits to shift to move 1 symbol
//#define BIT_PER_BYTE_MASK	0x0101010101010101

const int8_t gf8_exp[] = {	// not a multiple of 2 so duplicate entries + offset needed for easy wraparound of negatives
	1,2,4,3,6,7,5,
	1,2,4,3,6,7,5
};

const int8_t gf8_log[] = {	// log_0 undefined so dummy 0xF8 included to simplify indexing
	0xF8,0,1,3,2,6,4,5
};

const uint32_t gf8_rs_G_polys[] = {
	      01,	//0 symbols (dummy for indexing)
	     011,	//1 symbol
	    0132,	//2 symbols
	   01753,	//3 symbols
	  014775,	//4 symbols
	 0122313,	//5 symbols
	01576342	//6 symbols (identical to entries 1 thru 7 of the exponent table)
};


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
	if (a ==0 || b == 0)
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
	r1 &= 006666666666;
	r0 ^= r1;
	r2 &= 004444444444;
	r0 ^= r2;

	return gf8_poly_reduce(r0, of);
}

//TODO: the resulting polynomial is longer than the inputs, might have to rethink the types if this exceeds 7 terms
// also might remove the size args since the loop continue should handle it fine and it's only a max array length of 7 for both

//Assumes that result can never be longer than 7 terms, and the shorter polynomial is in q
// currently assuming the second multiplier is no more than 3 terms, this is probably not enough
// may have to consider rearranging as interleaved high and low degree symbols if it results in faster ops
uint32_t gf8_poly_mul(uint32_t p, uint32_t q)
{
	//TODO: check if splitting the terms in half makes the recombining faster
	uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, of;
	r0 = (q & 1) ? p : 0;
	r1 = (q & 2) * p;
	r2 = (q & 4) * p;
	r3 = (q & 8) * p;
	r4 = (q & 16) * p;
	r5 = (q & 32) * p;
	r6 = (q & 64) * p;
	r7 = (q & 128) * p;
	r8 = (q & 256) * p;
	of = (r1 & 011111111110) ^ (r2 & 033333333330);
	of ^= (r3 & 011111111110) ^ (r4 & 033333333330);
	of ^= (r5 & 011111111110) ^ (r6 & 033333333330);
	of ^= (r7 & 011111111110) ^ (r8 & 033333333330);
	r1 &= 006666666666;
	r0 ^= r1;
	r2 &= 004444444444;
	r0 ^= r2;
	r3 &= 006666666666;
	r0 ^= r3;
	r4 &= 004444444444;
	r0 ^= r4;
	r5 &= 006666666666;
	r0 ^= r5;
	r6 &= 004444444444;
	r0 ^= r6;
	r7 &= 006666666666;
	r0 ^= r7;
	r8 &= 004444444444;
	r0 ^= r8;

	return gf8_poly_reduce(r0, of);
}

//p is dividend, q is divisor, p_sz and q_sz are size in BITS not symbols
//returns remainder of the division since the quotient is never used
uint32_t gf8_poly_mod(uint32_t p, int8_t p_sz, uint32_t q, int8_t q_sz)
{
	//if p_sz and q_sz is known at compile time, this can be rewritten to be unrollable
	p_sz -= GF8_IDX_INC;
	q_sz -= GF8_IDX_INC;
	//uncomment the following line to return the quotient and remainder in a single Poly8 with the start of the quotient at b_arr[q_sz - 2]
	//q &= ~(-1 << q_sz); //clears the highest order term which should be a 1
	p = p << (8*q_sz);
	q = q << (8*p_sz);
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
    for(--p_sz; p_sz >= 0; p_sz -= GF8_IDX_INC)
        y = gf8_exp[gf8_log[y] + logx] ^ ((p >> p_sz) & 7);

    return y;
}
/*
//returns a generator polynomial for BCH view Reed Solomon with a given number of check symbols
uint32_t gf8_rs_G_unpack(int8_t syms)
{
	// packed octal format of the generator polynomials for given numbers of check symbols,
	// omits the trailing highest order 1 term 
	// indexing is done by shifting 3*n*(n-1)/2 where n is the number of check symbols
	// and then shifting and masking 3 bits at a time for a count of n
	uint64_t packed = 0576342223134775753321;
	/*
	polynomial coefficients for 1 thru 6 symbols highest to lowest term order are as follows
	          1,1	1 symbol
	        1,3,2	2 symbols
	      1,7,5,3	3 symbols
	    1,4,7,7,5	4 symbols
	  1,2,2,3,1,3	5 symbols
	1,5,7,6,3,4,2	6 symbols (identical to entries 1 thru 7 of the exponent table)
	
	int8_t chk_bits = syms * GF8_IDX_INC;
	packed >>= chk_bits * (syms-1) / 2;
	packed &= ~(-1 << chk_bits);
	return packed | (1 << chk_bits);
}
*/
//encodes a block of up to 18 bits worth of raw data as a Reed Solomon code word
//infers message length from provided data, doesn't verify that data length fits with the intended number of check symbols
//TODO: might be better to have the check polynomials easily indexable by check symbol count, consider rewriting for that
uint32_t gf8_rs_encode(int32_t raw, uint32_t chk_poly, int8_t chk_sz)
{
	int8_t msg_sz = 0;
	for(int8_t i = 1; i < raw; i << GF8_IDX_INC)
		msg_sz += GF8_IDX_INC;
	
	uint32_t chk = gf8_poly_mod(raw, msg_sz, chk_poly, chk_sz);
	raw <<= (chk_sz - GF8_IDX_INC);
	return raw | chk;
}

/*	//maybe do this later but GF8 is small enough that I did it by hand and can include all of them in a small space
//returns a generator polynomial for BCH view Reed Solomon with a given number of check symbols
uint32_t gf8_rs_gen_poly(int8_t chk_syms) {}
*/
/*
int8_t cl_div3_noLUT(int8_t x)
{
	while (x & 0xF8)
	{
		x ^= PRIME_GF8 << (28 - __builtin_clz(x));	//max of 3 bits shifted + 24 since builtin operates on uint32, not uint8
	}
	return x;
}

int8_t cl_mul3_noLUT(int8_t a, int8_t b)
{
	int8_t accum0, accum1, accum2;
	accum0 = (b & 1) ? a : 0;
	accum1 = (b & 2) ? a << 1 : 0;
	accum2 = (b & 4) ? a << 2 : 0;
	accum0 ^= accum1 ^ accum2;
	//printf("acc: %i\n", accum0);
	return 	cl_div3_noLUT(accum0);
}
*/