#ifndef RS_GF8_H
#define RS_GF8_H

#include <stdint.h>

extern const uint32_t rs8_G_polys[7];
//used to simplify memory allocation and vectorize some intstructions,
// only works for GF(8) because symbol length is limited to 7
/*uint32_t {
    uint64_t ull;
    int8_t b_arr[8];
};*/

uint32_t rs8_encode(uint32_t raw, uint32_t chk_poly, int8_t chk_sz);

uint32_t rs8_get_syndromes(uint32_t p, int8_t p_sz, int8_t nsyms);

uint32_t rs8_get_errata_locator(int8_t e_pos);

uint32_t rs8_get_errata_evaluator(uint32_t synd, int8_t chk_sz, uint32_t e_loc);
//uint32_t rs8_get_forney_syndromes(uint32_t synd, int8_t s_sz, uint32_t e_loc, int8_t e_sz);

uint32_t rs8_get_errata_magnitude(uint32_t synd, int8_t chk_sz, uint32_t e_loc, int8_t e_pos);

#endif  //RS_GF8_H