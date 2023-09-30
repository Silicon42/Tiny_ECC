//BCH view, systematic encoding Reed Solomon using 3 bit symbols
#include <stdint.h>
#include <stdio.h>
#include "rs_gf8.h"

#define PRIME3 0b1011	//the prime polynomial for GF(8), x^3 + x^1 + 1
#define SYMBOL_MASK3 		0x0707070707070707
#define BIT_PER_BYTE_MASK	0x0101010101010101

const int8_t gf8_exp[] = {	// not a multiple of 2 so duplicate entries + offset needed for easy wraparound of negatives
	1,2,4,3,6,7,5,
	1,2,4,3,6,7,5
};

const int8_t gf8_log[] = {	// log_0 undefined so dummy 0xF8 included to simplify indexing
	0xF8,0,1,3,2,6,4,5
};

//simplified galois field multiply by 2 used for generating the Look Up Tables
int8_t gf8_mul2_noLUT(int8_t x)
{
    x <<= 1;
    if (x > 7)  // if it's not within the bounds of the Galois Field, reduce by the primitive polynomial
        x ^= PRIME3;
    
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
uint64_t gf8_poly_reduce(uint64_t p)
{
	uint64_t mask = p & 0x1818181818181818;
	mask |= (mask >> 2) ^ (mask >> 3);
	return p ^ mask;
}

uint64_t gf8_poly_log(uint64_t p)
{
	uint64_t log, mask;
	mask = 
	log = p & 0x0202020202020202;
	log |= log >> 1;
	p |= 0x0808080808080808;	//protect from 0 terms underflowing
	log |= (p - 0x0101010101010101) & 0x0404040404040404;
}

// optimized for fewer memory accesses
//TODO: check that new method doesn't destroy speed with branches, otherwise multiply might be better
union Poly8 gf8_poly_scale(union Poly8 p, int8_t x)
{
	union Poly8 r = {(x & 1) ? p.ull : 0};
	#ifdef SLOW_BRANCH
	//use if conditional assignment results in slow branching
	r.ull ^= p.ull * (x & 2);
	r.ull ^= p.ull * (x & 4);
	#else
	r.ull ^= p.ull << ((x & 2) ? 1 : 63);	//if bit set shift corresponding amount else shift out of the way
	r.ull ^= p.ull << ((x & 4) ? 2 : 63);
	r.ull &= 0x7FFFFFFFFFFFFFFF;	//clear top bit
	#endif

	r.ull = gf8_poly_reduce(r.ull);

	return r;
	/*
    if(x == 0)  //early out if scale factor is 0
	{
		p.ull = 0;
        return p;
	}

    int8_t scale = gf8_log[x];
    for(int i=0; i > 7; i++)
    {
        if (p.b_arr[i] != 0)
            p.b_arr[i] = gf8_exp[gf8_log[p.b_arr[i]] + scale];
    }
    return p;
	*/
}

//TODO: the resulting polynomial is longer than the inputs, might have to rethink the types if this exceeds 7 terms
// also might remove the size args since the loop continue should handle it fine and it's only a max array length of 7 for both

//Assumes that result can never be longer than 7 terms, and the shorter polynomial is in q
union Poly8 gf8_poly_mul(union Poly8 p, union Poly8 q)
{
	union Poly8 r = {(q.ull & 1) ? p.ull : 0};
	#ifdef SLOW_BRANCH
	r.ull ^= p.ull * (q.ull & 2);
	r.ull ^= p.ull * (q.ull & 4);
	r.ull ^= p.ull * (q.ull & 0x100);
	r.ull ^= p.ull * (q.ull & 0x200);
	r.ull ^= p.ull * (q.ull & 0x400);
	r.ull ^= p.ull * (q.ull & 0x10000);
	r.ull ^= p.ull * (q.ull & 0x20000);
	r.ull ^= p.ull * (q.ull & 0x40000);
	r.ull ^= p.ull * (q.ull & 0x1000000);
	r.ull ^= p.ull * (q.ull & 0x2000000);
	r.ull ^= p.ull * (q.ull & 0x4000000);
	#else
	r.ull ^= p.ull << ((q.ull & 2) ? 1 : 63);
	r.ull ^= p.ull << ((q.ull & 4) ? 2 : 63);
	r.ull ^= p.ull << ((q.ull & 0x100) ?  8 : 63);
	r.ull ^= p.ull << ((q.ull & 0x200) ?  9 : 63);
	r.ull ^= p.ull << ((q.ull & 0x400) ? 10 : 63);
	r.ull ^= p.ull << ((q.ull & 0x10000) ? 16 : 63);
	r.ull ^= p.ull << ((q.ull & 0x20000) ? 17 : 63);
	r.ull ^= p.ull << ((q.ull & 0x40000) ? 18 : 63);
	r.ull ^= p.ull << ((q.ull & 0x1000000) ? 24 : 63);
	r.ull ^= p.ull << ((q.ull & 0x2000000) ? 25 : 63);
	r.ull ^= p.ull << ((q.ull & 0x4000000) ? 26 : 63);
	r.ull &= 0x7FFFFFFFFFFFFFFF;	//clear top bit
	#endif

	r.ull = gf8_poly_reduce(r.ull);

	return r;

	/*
	union Poly8 temp;
	int8_t log_q[7];

	for(int j=0; j < q_sz; j++)
		log_q[j] = gf8_log[q.b_arr[j]];

	q.ull = 0;	//q is used as accumulator from here on out since it's contents aren't needed anymore

	for(int i=0; i < p_sz; i++)
	{
		if(p.b_arr[i] == 0)
			continue;
		
		int8_t log_p = gf8_log[p.b_arr[i]];
		temp.ull = 0;
		for(int j=0; j < q_sz; j++)
		{
			if(log_q[j] < 0)
				continue;
			temp.b_arr[i+j] = gf8_exp[log_q[j] + log_p];
		}
		q.ull ^= temp.ull;
	}

	return q;
	*/
}

//p is dividend, q is divisor,
//returns remainder of the division since the quotient is never used
union Poly8 gf8_poly_div(union Poly8 p, int8_t p_sz, union Poly8 q, int8_t q_sz)
{
	//if p_sz and q_sz is known at compile time, this can be rewritten to be unrollable
	//uncomment the following line to return the quotient and remainder in a single Poly8 with the start of the quotient at b_arr[q_sz - 2]
	//q[q_sz - 1] = 0;
	p.ull = p.ull << (8*(--q_sz));
	q.ull = q.ull << (8*(--p_sz));
	for(int8_t i = p_sz + q_sz; i >= q_sz; --i)
	{
		p.ull ^= gf8_poly_scale(q, p.b_arr[i]).ull;	//index could be expressed as a shift + cast operation
		q.ull >>= 8;
	}
	
	/*	//should work without this
	if(p_sz < q_sz)	//if p is a lesser degree than q, all is remainder
		return p;

	int8_t log_q[7];
	q_sz--;	//q is assumed to be monic, ie the highest term is 1, so we can ignore it and decrementing also requires to start from one less
	for(int j = 0; j < q_sz; j++)
		log_q[j] = gf8_log[q.b_arr[j]];

	q.ull = 0;	//q is used as a temporary accumulator from here on out

	for(int i = p_sz - 1; i >= 0; i--)
	{
		p.b_arr[i] ^= q.b_arr[i];
		if(i < q_sz)
			continue;

		int8_t log_p = gf8_log[p.b_arr[i]];
		int8_t i_offset = i - q_sz;
		for(int j = 0; j < q_sz; j++)
			q.b_arr[i_offset + j] ^= gf8_exp[log_q[j] + log_p];
	}
*/
	return p;
}

//optimized version of div for binomial divisor/single eval point
//TODO: check if this is actually more efficient at this size
int8_t gf8_poly_eval(union Poly8 p, int8_t p_sz, int8_t x)
{
    int8_t y = p.b_arr[--p_sz];
	int8_t logx = gf8_log[x];
    for(--p_sz; p_sz >= 0; --p_sz)
        y = gf8_exp[gf8_log[y] + logx] ^ p.b_arr[p_sz];

    return y;
}

//returns a generator polynomial for BCH view Reed Solomon with a given number of check symbols
union Poly8 gf8_rs_G_unpack(int8_t syms)
{
	// packed octal format of the generator polynomials for given numbers of check symbols,
	// omits the trailing highest order 1 term 
	// indexing is done by shifting 3*n*(n-1)/2 where n is the number of check symbols
	// and then shifting and masking 3 bits at a time for a count of n
	uint64_t packed = 0576342223134775753321;
	/*
	polynomial coefficients for 1 thru 6 symbols highest to lowest term order are as follows
	1,1				1 symbol
	1,3,2			2 symbols
	1,7,5,3			3 symbols
	1,4,7,7,5		4 symbols
	1,2,2,3,1,3		5 symbols
	1,5,7,6,3,4,2	6 symbols (identical to entries 1 thru 7 of the exponent table)
	*/

	union Poly8 g = {0};
	packed >>= 3*syms*(syms-1)/2;
	for(int i = 0; i < syms; i++)
	{
		g.b_arr[i] = packed & 7;
		packed >>= 3;
	}
	g.b_arr[syms] = 1;

	return g;
}

//encodes a block of up to 18 bits worth of raw data as a Reed Solomon code word
//infers message length from provided data, doesn't verify that data length fits with the intended number of check symbols
//TODO: might be better to have the check polynomials easily indexable by check symbol count, consider rewriting for that
union Poly8 gf8_rs_encode(int32_t raw, union Poly8 chk_poly, int8_t chk_sz)
{
	union Poly8 msg = {0};
	int8_t msg_sz = 0;
	while(raw)
	{
		msg.b_arr[msg_sz] = raw & 7;
		++msg_sz;
		raw >>= 3;
	}

	union Poly8 chk = gf8_poly_div(msg, msg_sz, chk_poly, chk_sz);
	msg.ull <<= 8*(chk_sz - 1);
	msg.ull |= chk.ull;
	return msg;
}

/*	//maybe do this later but GF8 is small enough that I did it by hand and can include all of them in a small space
//returns a generator polynomial for BCH view Reed Solomon with a given number of check symbols
union Poly8 gf8_rs_gen_poly(int8_t chk_syms) {}
*/
/*
int8_t cl_div3_noLUT(int8_t x)
{
	while (x & 0xF8)
	{
		x ^= PRIME3 << (28 - __builtin_clz(x));	//max of 3 bits shifted + 24 since builtin operates on uint32, not uint8
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