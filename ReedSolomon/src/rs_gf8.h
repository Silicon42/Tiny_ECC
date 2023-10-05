#include "gf8.h"

//used to simplify memory allocation and vectorize some intstructions,
// only works for GF(8) because symbol length is limited to 7
/*uint32_t {
    uint64_t ull;
    int8_t b_arr[8];
};*/


uint32_t gf8_rs_G_unpack(int8_t syms);
