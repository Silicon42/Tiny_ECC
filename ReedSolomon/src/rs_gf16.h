#ifndef RS_GF16_H
#define RS_GF16_H

//BCH view, systematic encoding Reed Solomon using 4 bit symbols
//defined to use first consecutive root, c = 1 for slightly simplified decoding

#include <stdint.h>
#include "gf16.h"

gf16_poly rs16_encode(gf16_poly raw, int8_t chk_syms);

gf16_poly rs16_decode(gf16_poly recv, gf16_idx r_sz, int8_t chk_syms, int16_t e_pos, int16_t tx_pos);

#endif  //RS_GF16_H