#ifndef GF8_H
#define GF8_H

#include <stdint.h>

#define GF8_IDX_INC	3	    //how many bits to shift to move 1 symbol

//extern const int8_t gf8_exp[];	// length not a multiple of 2 so duplicate entries + offset needed for easy wraparound of negatives

extern const int8_t gf8_exp[14];
extern const int8_t gf8_log[8];
//extern const int8_t gf8_log[];	// log_0 undefined so dummy 0xF8 included to simplify indexing


int8_t gf8_mul2_noLUT(int8_t x);

int8_t gf8_mul(int8_t a, int8_t b);

int8_t gf8_div(int8_t a, int8_t b);

int8_t gf8_pow(int8_t x, int8_t power);

int8_t gf8_2pow(int8_t power);

int8_t gf8_inverse(int8_t x);

uint32_t gf8_poly_scale(uint32_t p, int8_t x);

uint32_t gf8_poly_mul(uint32_t p, uint32_t q);

int8_t gf8_poly_eval(uint32_t p, int8_t p_sz, int8_t x);

uint32_t gf8_poly_mod(uint32_t p, int8_t p_sz, uint32_t q, int8_t q_sz);

uint32_t gf8_poly_formal_derivative(uint32_t p);

#endif  //GF8_H