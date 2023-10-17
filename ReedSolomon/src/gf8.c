#include "gf8.h"
#include <stdio.h>  //TODO: remove this and the associated printf

#define PRIME_GF8 0b1011	//the prime polynomial for GF(8), x^3 + x^1 + 1
//masks that isolate out the term overflow from the result in the mul and scale functions
#define GF8_R1_OF 011111111110
#define GF8_R2_OF 033333333330
#define GF8_R1_R0 006666666666
#define GF8_R2_R0 004444444444
//mask to isolate just the odd terms for the formal derivative
#define GF8_ODD   007070707070

const int8_t gf8_exp[14] = {	// length not a multiple of 2 so duplicate entries + offset needed for easy wraparound of negatives
	1,2,4,3,6,7,5,
	1,2,4,3,6,7,5
};

const int8_t gf8_log[8] = {	// log_0 undefined so dummy 0xF8 included to simplify indexing
	0xF8,0,1,3,2,6,4,5
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
	if (a == 0 || b == 0)
		return 0;
	
	return gf8_exp[gf8_log[a] + gf8_log[b]];
}

int8_t gf8_pow(int8_t x, int8_t power)
{
	return gf8_exp[(gf8_log[x] * power) % 14];	// 14 is gf8_exp table size 
}

// slight optimization since most calls use x = 2 which evaluates to 1
// power is assumed to be in the range of 0 to 13
int8_t gf8_2pow(int8_t power)
{
	return gf8_exp[power];
}

int8_t gf8_inverse(int8_t x)
{
	return gf8_exp[7 - gf8_log[x]];
}

//prior to reduction, term can extend up to 2 bits above symbol due to shifting
//this function is customized to GF(8) with prime polynomial 1011
uint32_t gf8_poly_reduce(uint32_t p, uint32_t of)
{
	return p ^ (of >> 2) ^ (of >> 3);
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

	of  = (r1 & GF8_R1_OF) ^ (r2 & GF8_R2_OF);
	r0 ^= (r1 & GF8_R1_R0) ^ (r2 & GF8_R2_R0);

	return gf8_poly_reduce(r0, of);
}

//FIXME: the resulting polynomial is longer than the inputs and in Reed-Solomon with 6 check symbols can overflow
// if the result would be more than 10 terms long, need to add an alternate version that doesn't have this limitation
// for completness' sake of the Galois Field arithmetic

//Assumes that result can never be longer than 10 terms, and the shorter polynomial is in q
// currently assuming the second multiplier is no more than 5 terms, this is probably not enough
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

	r0 ^= (q & 01000) * p;
	r1 ^= (q & 02000) * p;
	r2 ^= (q & 04000) * p;

	r0 ^= (q & 010000) * p;
	r1 ^= (q & 020000) * p;
	r2 ^= (q & 040000) * p;

	r0 ^= (q & 0100000) * p;
	r1 ^= (q & 0200000) * p;
	r2 ^= (q & 0400000) * p;

	of  = (r1 & GF8_R1_OF) ^ (r2 & GF8_R2_OF);
	r0 ^= (r1 & GF8_R1_R0) ^ (r2 & GF8_R2_R0);

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
	//q &= ~((uint32_t)-1 << q_sz); //clears the highest order term which should be a 1
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
    {
        if(y)
            y = gf8_exp[gf8_log[y] + logx];
        
        y ^= ((p >> p_sz) & 7);
    }
    return y;
}

//formal derivative of characteristic 2 keeps only the odd polynomials and reduces the degree by 1 step
uint32_t gf8_poly_formal_derivative(uint32_t p)
{
    return (p & GF8_ODD) >> GF8_IDX_INC;
}