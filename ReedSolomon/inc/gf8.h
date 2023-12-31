#ifndef GF8_H
#define GF8_H

// GF(8) field element and polynomial math library
//
// polynomials are packed as consecutive octals in uint32 with lowest order terms in the least significant bits
//  as such they are limited to a max of 10 terms although polynomial multiply is further limited to only
//  functioning correctly for up to 6 terms in the second argument
//
// GF(8) as defined by using primitive polynomial 1x^3 + 0x^2 + 1x^1 + 1x^0, ie 1011 for field element reduction
//  and the primitive field element used is 2 for the exponent and log tables

#include <stdint.h>

#define GF8_SYM_SZ 3				// how many bits to shift to move 1 symbol
#define GF8_MAX 7					// max value a field element can have
#define GF8_EXP_ENTRIES 2 * GF8_MAX // number of entries in the exponent table

typedef int8_t gf8_idx;		// represents a polynomial term index or size in terms of bits, should always be incremented/decremented by GF8_SYM_SZ
typedef int8_t gf8_elem;	// a single GF(8) element, only valid in the range of 0 through 7
typedef int32_t gf8_poly;	// GF(8) polynomial of order no greater than 9 (10 terms) packed in a uint32

extern const gf8_elem gf8_exp[GF8_EXP_ENTRIES];	// length not a multiple of 2 so duplicate entries + offset needed for fast wraparound of negatives
extern const gf8_elem gf8_log[1 + GF8_MAX];		// log_0 undefined so dummy 0xFF included to simplify indexing

gf8_elem gf8_mul2_noLUT(gf8_elem x);

gf8_elem gf8_mul(gf8_elem a, gf8_elem b);

gf8_elem gf8_div(gf8_elem a, gf8_elem b);

gf8_elem gf8_pow(gf8_elem x, int8_t power);

gf8_elem gf8_2pow(int8_t power);

gf8_elem gf8_inverse(gf8_elem x);

gf8_poly gf8_poly_scale(gf8_poly p, gf8_elem x);

gf8_poly gf8_poly_mul(gf8_poly p, gf8_poly q);

gf8_poly gf8_poly_mul_q0_monic(gf8_poly p, gf8_poly q);

gf8_elem gf8_poly_eval(gf8_poly p, gf8_idx p_sz, gf8_elem x);

gf8_poly gf8_poly_mod(gf8_poly p, gf8_idx p_sz, gf8_poly q, gf8_idx q_sz);

gf8_poly gf8_poly_formal_derivative(gf8_poly p);

int8_t gf8_poly_get_order(gf8_poly p);

gf8_idx gf8_poly_get_size(gf8_poly p);

#endif // GF8_H