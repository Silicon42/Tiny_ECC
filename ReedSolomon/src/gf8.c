#include "gf8.h"

#define PRIME_GF8 0b1011 // the prime polynomial for GF(8), x^3 + x^1 + 1
// masks that isolate out the term overflow from the result in the mul and scale functions
#define GF8_R1_OF 011111111110
#define GF8_R2_OF 033333333330
#define GF8_R1_R0 006666666666
#define GF8_R2_R0 004444444444
// mask to isolate just the odd terms for the formal derivative
#define GF8_ODD   007070707070

const gf8_elem gf8_exp[GF8_EXP_ENTRIES] = {	// length not a multiple of 2 so duplicate entries + offset needed for easy wraparound of negatives
	1, 2, 4, 3, 6, 7, 5,
	1, 2, 4, 3, 6, 7, 5};

const gf8_elem gf8_log[8] = {	// log_0 undefined so dummy 0xFF included to simplify indexing
	0xFF, 0, 1, 3, 2, 6, 4, 5};

// simplified galois field multiply by 2 used for generating the Look Up Tables
gf8_elem gf8_mul2_noLUT(gf8_elem x)
{
	x <<= 1;
	if (x > GF8_MAX)	// if it's not within the bounds of the Galois Field, reduce by the primitive polynomial
		x ^= PRIME_GF8;

	return x;
}

gf8_elem gf8_div(gf8_elem a, gf8_elem b)
{
	if (b == 0)
		return -1;	// divide by 0 error, normal operation should never get here
	if (a == 0)
		return 0;

	return gf8_exp[gf8_log[a] - gf8_log[b] + GF8_MAX];	// +GF8_MAX offset to keep range positive
}

gf8_elem gf8_mul(gf8_elem a, gf8_elem b)
{
	if (a == 0 || b == 0)
		return 0;

	return gf8_exp[gf8_log[a] + gf8_log[b]];
}

gf8_elem gf8_pow(gf8_elem x, int8_t power)
{
	return gf8_exp[(gf8_log[x] * power) % GF8_EXP_ENTRIES];
}

// slight optimization since most calls use x = 2 which evaluates to 1
// power is assumed to be in the range of 0 to 13
gf8_elem gf8_2pow(int8_t power)
{
	return gf8_exp[power];
}

gf8_elem gf8_inverse(gf8_elem x)
{
	return gf8_exp[GF8_MAX - gf8_log[x]];
}

// prior to reduction, term can extend up to 2 bits above symbol due to shifting
// this function is customized to GF(8) with prime polynomial 1011
gf8_poly gf8_poly_reduce(gf8_poly p, gf8_poly of)
{
	return p ^ (of >> 2) ^ (of >> 3);
}

// optimized for fewer memory accesses
// TODO: check if multiplies are faster, currently assuming that single shifts and conditional assignment are better
gf8_poly gf8_poly_scale(gf8_poly p, gf8_elem x)
{
	gf8_poly r0, r1, r2, of;
	r0 = (x & 1) ? p : 0;
	p <<= 1;
	r1 = (x & 2) ? p : 0;
	p <<= 1;
	r2 = (x & 4) ? p : 0;

	of = (r1 & GF8_R1_OF) ^ (r2 & GF8_R2_OF);
	r0 ^= (r1 & GF8_R1_R0) ^ (r2 & GF8_R2_R0);

	return gf8_poly_reduce(r0, of);
}

// Assumes that result can never be longer than 10 terms, and the shorter polynomial is in q
//  currently assuming the second multiplier is no more than 5 terms, this is just enough for
//  Reed Solomon with a max of 6 check symbols with specific optimizations
gf8_poly gf8_poly_mul(gf8_poly p, gf8_poly q)
{
	gf8_poly r0, r1, r2, of;
	// term 0
	r0 = (q & 01) ? p : 0;
	r1 = (q & 02) * p;
	r2 = (q & 04) * p;
	// term 1
	r0 ^= (q & 010) * p;
	r1 ^= (q & 020) * p;
	r2 ^= (q & 040) * p;
	// term 2
	r0 ^= (q & 0100) * p;
	r1 ^= (q & 0200) * p;
	r2 ^= (q & 0400) * p;
	// term 3
	r0 ^= (q & 01000) * p;
	r1 ^= (q & 02000) * p;
	r2 ^= (q & 04000) * p;
	// term 4
	r0 ^= (q & 010000) * p;
	r1 ^= (q & 020000) * p;
	r2 ^= (q & 040000) * p;
	/*
		r0 ^= (q & 0100000) * p;
		r1 ^= (q & 0200000) * p;
		r2 ^= (q & 0400000) * p;
	*/
	of = (r1 & GF8_R1_OF) ^ (r2 & GF8_R2_OF);
	r0 ^= (r1 & GF8_R1_R0) ^ (r2 & GF8_R2_R0);

	return gf8_poly_reduce(r0, of);
}

// squeezes one more term out of poly_mul with the assumption that term 0 of q is always 1
gf8_poly gf8_poly_mul_q0_monic(gf8_poly p, gf8_poly q)
{
	return p ^ (gf8_poly_mul(p, q >> GF8_SYM_SZ) << GF8_SYM_SZ);
}

// p is dividend, q is divisor, p_sz and q_sz are size in BITS not symbols
// returns remainder of the division since the quotient is never used
gf8_poly gf8_poly_mod(gf8_poly p, gf8_idx p_sz, gf8_poly q, gf8_idx q_sz)
{
	// if p_sz and q_sz is known at compile time, this can be rewritten to be unrollable
	p_sz -= GF8_SYM_SZ;
	q_sz -= GF8_SYM_SZ;
	// uncomment the following line to return the quotient and remainder in a single return value with the start of the quotient at b_arr[q_sz - 2]
	// q &= ~((gf8_poly)-1 << q_sz); //clears the highest order term which should be a 1
	p <<= q_sz;
	q <<= p_sz;
	for (gf8_idx i = p_sz + q_sz; i >= q_sz; i -= GF8_SYM_SZ)
	{
		p ^= gf8_poly_scale(q, (p >> i) & GF8_MAX);
		q >>= GF8_SYM_SZ;
	}

	return p;
}

// optimized version of div for binomial divisor/single eval point
// TODO: check if this is actually more efficient at this size
gf8_elem gf8_poly_eval(gf8_poly p, gf8_idx p_sz, gf8_elem x)
{
	p_sz -= GF8_SYM_SZ;
	gf8_elem y = p >> p_sz;
	gf8_elem logx = gf8_log[x];
	for (p_sz -= GF8_SYM_SZ; p_sz >= 0; p_sz -= GF8_SYM_SZ)
	{
		if (y)
			y = gf8_exp[gf8_log[y] + logx];

		y ^= ((p >> p_sz) & GF8_MAX);
	}
	return y;
}

// formal derivative of characteristic 2 keeps only the odd polynomials and reduces the degree by 1 step
gf8_poly gf8_poly_formal_derivative(gf8_poly p)
{
	return (p & GF8_ODD) >> GF8_SYM_SZ;
}

int8_t gf8_poly_get_order(gf8_poly p)
{
	int8_t n = -1;
	for (gf8_poly i = 1; i <= p; i <<= GF8_SYM_SZ)
		++n;

	return n;
}

gf8_idx gf8_poly_get_size(gf8_poly p)
{
	gf8_idx p_sz = 0;
	for (gf8_poly i = 1; i <= p; i <<= GF8_SYM_SZ)
		p_sz += GF8_SYM_SZ;

	return p_sz;
}