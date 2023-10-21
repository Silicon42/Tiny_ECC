#ifndef RS_GF8_H
#define RS_GF8_H

//BCH view, systematic encoding Reed Solomon using 3 bit symbols
//defined to use first consecutive root = 1 for slightly simplified decoding

#include <stdint.h>
#include "gf8.h"

extern const uint32_t rs8_G_polys[7];
//used to simplify memory allocation and vectorize some intstructions,
// only works for GF(8) because symbol length is limited to 7
/*uint32_t {
    uint64_t ull;
    int8_t b_arr[8];
};*/

uint32_t rs8_encode(uint32_t raw, uint32_t chk_poly, gf8_idx chk_sz);

uint32_t rs8_decode(uint32_t recv, gf8_idx r_sz, int8_t chk_syms, int8_t e_pos);

uint32_t rs8_get_syndromes(uint32_t p, gf8_idx p_sz, int8_t nsyms);

uint32_t rs8_get_erasure_locator(int8_t e_pos);

uint32_t rs8_get_error_locator(uint32_t synd, gf8_idx s_sz);

uint32_t rs8_get_errata_evaluator(uint32_t synd, gf8_idx chk_sz, uint32_t e_loc);

int8_t rs8_get_error_pos(uint32_t error_loc, int8_t erase_pos);

uint32_t rs8_get_errata_magnitude(uint32_t e_eval, gf8_idx chk_sz, uint32_t e_loc, int8_t e_pos);

#endif  //RS_GF8_H