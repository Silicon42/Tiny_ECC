#include <stdint.h>

//used to simplify memory allocation and vectorize some intstructions,
// only works for GF(16) because symbol length is limited to 7
union Poly16 {
    __uint128_t uxl;
    uint64_t ull[2];
    int8_t b_arr[16];
};

int8_t gf16_mul2_noLUT(int8_t x);

int8_t gf16_mul(int8_t a, int8_t b);

int8_t gf16_div(int8_t a, int8_t b);

int8_t gf16_pow(int8_t x, int8_t power);

int8_t gf16_2pow(int8_t power);

int8_t gf16_inverse(int8_t x);

union Poly16 gf16_poly_scale(union Poly16 p, int8_t x);

union Poly16 gf16_poly_mul(union Poly16 p, union Poly16 q);

int8_t gf16_poly_eval(union Poly16 p, int8_t x);

union Poly16 gf16_poly_div(union Poly16 p, int8_t p_sz, union Poly16 q, int8_t q_sz);

union Poly16 gf16_rs_G_unpack(int8_t syms);
