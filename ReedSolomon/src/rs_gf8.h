#ifndef RS_GF8_H
#define RS_GF8_H

//BCH view, systematic encoding Reed Solomon using 3 bit symbols
//defined to use first consecutive root = 1 for slightly simplified decoding

#include <stdint.h>
#include "gf8.h"

gf8_poly rs8_encode(gf8_poly raw, int8_t chk_syms);

gf8_poly rs8_decode(gf8_poly recv, gf8_idx r_sz, int8_t chk_syms, int8_t e_pos, int8_t tx_pos);

#endif  //RS_GF8_H