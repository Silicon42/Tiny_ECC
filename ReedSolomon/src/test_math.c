#include "rs_gf8.h"
#include "gf8.h"
#include <stdio.h>

int main()
{
    uint32_t r;

    //TODO: make sure to test the individual symbol functions as well
    //TODO: do rigorous performance tests on various systems to get an idea which implementation of some functions
    // should be the default, ie multiplication vs shift + conditional assignement, for loop with term count tracking
    // vs term counting as needed vs max expected width mul/div to avoid branching, should there be an option for SIMD, etc.

    r = gf8_poly_scale(0137, 3);
    printf("%o\n", r);  //result 352

    r = gf8_poly_mul(0137, 035);
    printf("%o\n", r);  //result 3066

    r = gf8_poly_eval(017131, 15, 1);
    printf("%o\n", r);  //result 5

    r = gf8_poly_mod(01117, 12, 0132, 9);
    printf("%o\n", r);  //result 35

    r = rs8_encode(01117, 0132, 9);
    printf("%o\n", r);  //result 11735

    r = rs8_get_syndromes(0111700, 18, 2);
    printf("%o\n", r);  //result 36

    r = rs8_get_errata_locator(0b01000001);
    printf("%o\n", r);  //result 132

    return 0;
}