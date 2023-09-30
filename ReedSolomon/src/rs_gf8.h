#include <stdint.h>

//used to simplify memory allocation and vectorize some intstructions,
// only works for GF(8) because symbol length is limited to 7
union Poly8 {
    uint64_t ull;
    int8_t b_arr[8];
};

int8_t gf8_mul2_noLUT(int8_t x);

int8_t gf8_mul(int8_t a, int8_t b);

int8_t gf8_div(int8_t a, int8_t b);

int8_t gf8_pow(int8_t x, int8_t power);

int8_t gf8_2pow(int8_t power);

int8_t gf8_inverse(int8_t x);

union Poly8 gf8_poly_scale(union Poly8 p, int8_t x);

union Poly8 gf8_poly_mul(union Poly8 p, union Poly8 q);

int8_t gf8_poly_eval(union Poly8 p, int8_t p_sz, int8_t x);

union Poly8 gf8_poly_div(union Poly8 p, int8_t p_sz, union Poly8 q, int8_t q_sz);

union Poly8 gf8_rs_G_unpack(int8_t syms);
