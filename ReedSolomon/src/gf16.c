#include "gf16.h"

#define PRIME_GF16 0b10011	// the prime polynomial for GF(16), x^4 + x^1 + 1
// masks that isolate out the term overflow from the result in the mul and scale functions
#define GF16_R1_OF 0x1111111111111111
#define GF16_R2_OF 0x3333333333333333
#define GF16_R3_OF 0x7777777777777777
#define GF16_R1_R0 ~GF16_R1_OF
#define GF16_R2_R0 ~GF16_R2_OF
#define GF16_R3_R0 ~GF16_R3_OF
// mask to isolate just the odd terms for the formal derivative
#define GF16_ODD   0xF0F0F0F0F0F0F0F0

const gf16_elem gf16_exp[GF16_EXP_ENTRIES] = {	// length not a multiple of 2 so duplicate entries + offset needed for easy wraparound of negatives
	0x1, 0x2, 0x4, 0x8, 0x3, 0x6, 0xC, 0xB, 0x5, 0xA, 0x7, 0xE, 0xF, 0xD, 0x9,
	0x1, 0x2, 0x4, 0x8, 0x3, 0x6, 0xC, 0xB, 0x5, 0xA, 0x7, 0xE, 0xF, 0xD, 0x9};

const gf16_elem gf16_log[1 + GF16_MAX] = {	// log_0 undefined so dummy 0xFF included to simplify indexing
	0xFF, 0x0, 0x1, 0x4, 0x2, 0x8, 0x5, 0xA, 0x3, 0xE, 0x9, 0x7, 0x6, 0xD, 0xB, 0xC};

// simplified galois field multiply by 2 used for generating the Look Up Tables
gf16_elem gf16_mul2_noLUT(gf16_elem x)
{
	x <<= 1;
	if (x > GF16_MAX)	// if it's not within the bounds of the Galois Field, reduce by the primitive polynomial
		x ^= PRIME_GF16;

	return x;
}

gf16_elem gf16_div(gf16_elem a, gf16_elem b)
{
	if (b == 0)
		return -1;	// divide by 0 error, normal operation should never get here
	if (a == 0)
		return 0;

	return gf16_exp[gf16_log[a] - gf16_log[b] + GF16_MAX];	// +GF16_MAX offset to keep range positive
}

gf16_elem gf16_mul(gf16_elem a, gf16_elem b)
{
	if (a == 0 || b == 0)
		return 0;

	return gf16_exp[gf16_log[a] + gf16_log[b]];
}

gf16_elem gf16_pow(gf16_elem x, int8_t power)
{
	return gf16_exp[(gf16_log[x] * power) % GF16_EXP_ENTRIES];
}

// slight optimization since most calls use x = 2 which evaluates to 1
// power is assumed to be in the range of 0 to GF16_EXP_ENTRIES -1
gf16_elem gf16_2pow(int8_t power)
{
	return gf16_exp[power];
}

gf16_elem gf16_inverse(gf16_elem x)
{
	return gf16_exp[GF16_MAX - gf16_log[x]];
}

// prior to reduction, term can extend up to 2 bits above symbol due to shifting
// this function is customized to GF(16) with prime polynomial 10011
gf16_poly gf16_poly_reduce(gf16_poly p, gf16_poly of)
{
	return p ^ (of >> 3) ^ (of >> 4);
}

// optimized for fewer memory accesses
// TODO: check if multiplies are faster, currently assuming that single shifts and conditional assignment are better
gf16_poly gf16_poly_scale(gf16_poly p, gf16_elem x)
{
	gf16_poly r0, r1, r2, r3, of;
	r0 = (x & 1) ? p : 0;
	p <<= 1;
	r1 = (x & 2) ? p : 0;
	p <<= 1;
	r2 = (x & 4) ? p : 0;
	p <<= 1;
	r3 = (x & 8) ? p : 0;

	of = (r1 & GF16_R1_OF) ^ (r2 & GF16_R2_OF) ^ (r3 & GF16_R3_OF);
	r0 ^= (r1 & GF16_R1_R0) ^ (r2 & GF16_R2_R0) ^ (r3 & GF16_R3_R0);

	return gf16_poly_reduce(r0, of);
}

// FIXME: This is not enough terms for full GF(16) polynomials, it needs to be increased and
//  maybe converted to a loop
// Assumes that result can never be longer than 10 terms, and the shorter polynomial is in q
//  currently assuming the second multiplier is no more than 5 terms, this is just enough for
//  Reed Solomon with a max of 6 check symbols with specific optimizations
gf16_poly gf16_poly_mul(gf16_poly p, gf16_poly q)
{
	gf16_poly r0, r1, r2, r3, of;
	// term 0
	r0 = (q & 1) ? p : 0;
	r1 = (q & 2) * p;
	r2 = (q & 4) * p;
	r3 = (q & 8) * p;
	// term 1
	r0 ^= (q & 0x10) * p;
	r1 ^= (q & 0x20) * p;
	r2 ^= (q & 0x40) * p;
	r3 ^= (q & 0x80) * p;
	// term 2
	r0 ^= (q & 0x100) * p;
	r1 ^= (q & 0x200) * p;
	r2 ^= (q & 0x400) * p;
	r3 ^= (q & 0x800) * p;
	// term 3
	r0 ^= (q & 0x1000) * p;
	r1 ^= (q & 0x2000) * p;
	r2 ^= (q & 0x4000) * p;
	r3 ^= (q & 0x8000) * p;
	// term 4
	r0 ^= (q & 0x10000) * p;
	r1 ^= (q & 0x20000) * p;
	r2 ^= (q & 0x40000) * p;
	r3 ^= (q & 0x80000) * p;
	// term 5
	r0 ^= (q & 0x100000) * p;
	r1 ^= (q & 0x200000) * p;
	r2 ^= (q & 0x400000) * p;
	r3 ^= (q & 0x800000) * p;

	of = (r1 & GF16_R1_OF) ^ (r2 & GF16_R2_OF) ^ (r3 & GF16_R3_OF);
	r0 ^= (r1 & GF16_R1_R0) ^ (r2 & GF16_R2_R0) ^ (r3 & GF16_R3_R0);

	return gf16_poly_reduce(r0, of);
}

// squeezes one more term out of poly_mul with the assumption that term 0 of q is always 1
gf16_poly gf16_poly_mul_q0_monic(gf16_poly p, gf16_poly q)
{
	return p ^ (gf16_poly_mul(p, q >> GF16_SYM_SZ) << GF16_SYM_SZ);
}

// p is dividend, q is divisor, p_sz and q_sz are size in BITS not symbols
// returns remainder of the division since the quotient is never used
gf16_poly gf16_poly_mod(gf16_poly p, gf16_idx p_sz, gf16_poly q, gf16_idx q_sz)
{
	// if p_sz and q_sz is known at compile time, this can be rewritten to be unrollable
	p_sz -= GF16_SYM_SZ;
	q_sz -= GF16_SYM_SZ;
	// uncomment the following line to return the quotient and remainder in a single return value with the start of the quotient at b_arr[q_sz - 2]
	// q &= ~((gf16_poly)-1 << q_sz); //clears the highest order term which should be a 1
	p <<= q_sz;
	q <<= p_sz;
	for (gf16_idx i = p_sz + q_sz; i >= q_sz; i -= GF16_SYM_SZ)
	{
		p ^= gf16_poly_scale(q, (p >> i) & GF16_MAX);
		q >>= GF16_SYM_SZ;
	}

	return p;
}

// optimized version of div for binomial divisor/single eval point
// TODO: check if this is actually more efficient at this size
gf16_elem gf16_poly_eval(gf16_poly p, gf16_idx p_sz, gf16_elem x)
{
	p_sz -= GF16_SYM_SZ;
	gf16_elem y = p >> p_sz;
	gf16_elem logx = gf16_log[x];
	for (p_sz -= GF16_SYM_SZ; p_sz >= 0; p_sz -= GF16_SYM_SZ)
	{
		if (y)
			y = gf16_exp[gf16_log[y] + logx];

		y ^= ((p >> p_sz) & GF16_MAX);
	}
	return y;
}

// formal derivative of characteristic 2 keeps only the odd polynomials and reduces the degree by 1 step
gf16_poly gf16_poly_formal_derivative(gf16_poly p)
{
	return (p & GF16_ODD) >> GF16_SYM_SZ;
}

int8_t gf16_poly_get_order(gf16_poly p)
{
	int8_t n = -1;
	for (gf16_poly i = 1; i <= p; i <<= GF16_SYM_SZ)
		++n;

	return n;
}

gf16_idx gf16_poly_get_size(gf16_poly p)
{
	gf16_idx p_sz = 0;
	for (gf16_poly i = 1; i <= p; i <<= GF16_SYM_SZ)
		p_sz += GF16_SYM_SZ;

	return p_sz;
}