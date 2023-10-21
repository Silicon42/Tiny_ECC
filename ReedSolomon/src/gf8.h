#ifndef GF8_H
#define GF8_H

//GF(8) field element and polynomial math library
//
//polynomials are packed as consecutive octals in uint32 with lowest order terms in the least significant bits
// as such they are limited to a max of 10 terms although polynomial multiply is further limited to only
// functioning correctly for upt to 6 terms in the second argument
//
//GF(8) as defined by using primitive polynomial 1x^3 + 0x^2 + 1x^1 + 1x^0, ie 1011 for field element reduction
// and the primitive field element used is 2 for the exponent and log tables


#include <stdint.h>

#define GF8_SYM_SZ  3	    //how many bits to shift to move 1 symbol

extern const int8_t gf8_exp[14];	//length not a multiple of 2 so duplicate entries + offset needed for fast wraparound of negatives
extern const int8_t gf8_log[8];	    //log_0 undefined so dummy 0xF8 included to simplify indexing

typedef int8_t gf8_idx;     //represents a polynomial term index or size in terms of bits, should always be incremented/decremented via functions, never operators
typedef int8_t gf8_elem;    //a single GF(8) element
typedef uint32_t gf8_poly;  //a packed GF(8) polynomial of order no greater than 9 (10 terms)
//FIXME: finish converting to new types for clarity and type safety

int8_t gf8_idx_inc(gf8_idx i);
int8_t gf8_idx_dec(gf8_idx i);

//inline gf8_idx gf8_idx_inc(gf8_idx i) {return i + GF8_SYM_SZ;}
//inline gf8_idx gf8_idx_dec(gf8_idx i) {return i - GF8_SYM_SZ;}

int8_t gf8_mul2_noLUT(int8_t x);

int8_t gf8_mul(int8_t a, int8_t b);

int8_t gf8_div(int8_t a, int8_t b);

int8_t gf8_pow(int8_t x, int8_t power);

int8_t gf8_2pow(int8_t power);

int8_t gf8_inverse(int8_t x);

uint32_t gf8_poly_scale(uint32_t p, int8_t x);

uint32_t gf8_poly_mul(uint32_t p, uint32_t q);

int8_t gf8_poly_eval(uint32_t p, gf8_idx p_sz, int8_t x);

uint32_t gf8_poly_mod(uint32_t p, gf8_idx p_sz, uint32_t q, gf8_idx q_sz);

uint32_t gf8_poly_formal_derivative(uint32_t p);

int8_t gf8_poly_get_order(uint32_t p);

#endif  //GF8_H