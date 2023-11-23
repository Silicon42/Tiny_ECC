#ifndef GF16_H
#define GF16_H

// GF(16) field element and polynomial math library
//
// polynomials are packed as consecutive octals in uint64 with lowest order terms in the least significant bits
//  as such they are limited to a max of 10 terms although polynomial multiply is further limited to only
//  functioning correctly for up to 6 terms in the second argument
//
// GF(16) as defined by using primitive polynomial 1x^3 + 0x^2 + 1x^1 + 1x^0, ie 1011 for field element reduction
//  and the primitive field element used is 2 for the exponent and log tables

#include <stdint.h>

#define GF16_SYM_SZ 4					// how many bits to shift to move 1 symbol
#define GF16_MAX 15						// max value a field element can have
#define GF16_EXP_ENTRIES 2 * GF16_MAX	// number of entries in the exponent table

typedef int8_t gf16_idx;	// represents a polynomial term index or size in terms of bits, should always be incremented/decremented by GF16_SYM_SZ
typedef int8_t gf16_elem;	// a single GF(16) element, only valid in the range of 0 through 7
typedef int64_t gf16_poly;	// GF(16) polynomial of order no greater than 14 (15 terms) packed in a uint64,
// while there is room for a 16th term, there is not for its overflow and BCH view Reed Solomon is limited to 15 anyway

extern const gf16_elem gf16_exp[GF16_EXP_ENTRIES];	// length not a multiple of 2 so duplicate entries + offset needed for fast wraparound of negatives
extern const gf16_elem gf16_log[1 + GF16_MAX];		// log_0 undefined so dummy 0xFF included to simplify indexing

gf16_elem gf16_mul2_noLUT(gf16_elem x);

gf16_elem gf16_mul(gf16_elem a, gf16_elem b);

gf16_elem gf16_div(gf16_elem a, gf16_elem b);

gf16_elem gf16_pow(gf16_elem x, int8_t power);

gf16_elem gf16_2pow(int8_t power);

gf16_elem gf16_inverse(gf16_elem x);

gf16_poly gf16_poly_scale(gf16_poly p, gf16_elem x);

gf16_poly gf16_poly_mul(gf16_poly p, gf16_poly q);

gf16_poly gf16_poly_mul_q0_monic(gf16_poly p, gf16_poly q);

gf16_elem gf16_poly_eval(gf16_poly p, gf16_idx p_sz, gf16_elem x);

gf16_poly gf16_poly_mod(gf16_poly p, gf16_idx p_sz, gf16_poly q, gf16_idx q_sz);

gf16_poly gf16_poly_formal_derivative(gf16_poly p);

int8_t gf16_poly_get_order(gf16_poly p);

gf16_idx gf16_poly_get_size(gf16_poly p);

#endif // GF16_H