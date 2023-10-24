#include "rs_gf8.h"
#include "gf8.h"
#include <stdio.h>

int main()
{
    gf8_poly r, synd, e_loc, e_eval;

    //TODO: make sure to test the individual symbol functions as well
    //TODO: do rigorous performance tests on various systems to get an idea which implementation of some functions
    // should be the default, ie multiplication vs shift + conditional assignement, for loop with term count tracking
    // vs term counting as needed vs max expected width mul/div to avoid branching, should there be an option for SIMD, etc.
/*
    r = gf8_poly_scale(0137, 3);
    printf("%o\n", r);  //result 352

    r = gf8_poly_mul(0137, 035);
    printf("%o\n", r);  //result 3066

    r = gf8_poly_eval(017131, 15, 1);
    printf("%o\n", r);  //result 5

    r = gf8_poly_mod(01117, 12, 0132, 9);
    printf("%o\n", r);  //result 35

    r = gf8_poly_mod(01117, 15, 0132, 9);   //check that longer than true p_sz works
    printf("%o\n", r);  //result 35

    printf("\n");
    //2 check symbols, 4 data symbols
    r = rs8_encode(01117, 0132, 9);
    printf("%o\n", r);  //result 11735

    printf("\n");
    //2 check symbols, 2 erasures
    synd = rs8_get_syndromes(0111700, 18, 2);
    printf("%o\n", synd);  //result 36

    e_loc = rs8_get_errata_locator(0b00000011);
    printf("%o\n", e_loc);  //result 145

    //uint32_t m_synd = rs8_get_forney_syndromes(synd, 6, e_loc, 9);
    //printf("%o\n", m_synd);  //result 36

    printf("\n");
    //2 check symbols, 1 erasure
    synd = rs8_get_syndromes(0111705, 18, 2);
    printf("%o\n", synd);  //result 63

    e_loc = rs8_get_errata_locator(0b00000010);
    printf("%o\n", e_loc);  //result 15

    //m_synd = rs8_get_forney_syndromes(synd, 6, e_loc, 6);
    //printf("%o\n", m_synd);  //result 60

    printf("\n");
    //6 check symbols, 1 data symbol
    r = rs8_encode(01, 01576342, 21);
    printf("%o\n", r);  //result 1576342

    //6 check symbols, 4 erasures
    synd = rs8_get_syndromes(0342, 21, 6);  //also tests that longer than true p_sz works
    printf("%o\n", synd);  //result 611565

    e_loc = rs8_get_errata_locator(0b01111000);
    printf("%o\n", e_loc);  //result 13123

    //m_synd = rs8_get_forney_syndromes(synd, 21, e_loc, 18);
    //printf("%o\n", m_synd);  //result 611500

    r = rs8_get_errata_magnitude(synd, 21, e_loc, 0b01111000);
    printf("%o\n", r);  //result 1576000
*/
    //4 check symbols, 3 data symbols
    r = rs8_encode(0123, rs8_G_polys[4], 15);
    printf("%o\n", r);  //result 1230013

    //4 check symbols, 2 errors
    r = rs8_decode(030013, 21, 4, 0, 0x7F);
    printf("%o\n", r);  //result 1230013

    //4 check symbols, 2 erasures
    r = rs8_decode(030013, 21, 4, 0b1100000, 0x7F);
    printf("%o\n", r);  //result 1230013

    //4 check symbols, 3 erasures
    r = rs8_decode(00013, 21, 4, 0b1110000, 0x7F);
    printf("%o\n", r);  //result 1230013

    //4 check symbols, 4 erasures
    r = rs8_decode(01013, 21, 4, 0b1111000, 0x7F);
    printf("%o\n", r);  //result 1230013

    //4 check symbols, 2 erasures, 1 error
    r = rs8_decode(00013, 21, 4, 0b1100000, 0x7F);
    printf("%o\n", r);  //result 1230013

    //4 check symbols, no errata
    r = rs8_decode(01230013, 21, 4, 0, 0x7F);
    printf("%o\n", r);  //result 1230013
    
    return 0;
}