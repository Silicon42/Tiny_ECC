#include <stdint.h>

//used to simplify memory allocation and vectorize some intstructions,
// only works for GF(8) because symbol length is limited to 7
/*uint32_t {
    uint64_t ull;
    int8_t b_arr[8];
};*/

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

uint32_t gf8_rs_G_unpack(int8_t syms);
