#include "rs_gf8.h"
#include "gf8.h"
#include <stdio.h>

int main()
{
    uint32_t p, q, r;
    p = 0137;
    q = 035;

    r = gf8_poly_mul(p, q);

    printf("%o\n", r);

   // r = gf8_rs_G_unpack(3);

    //printf("\n%llX", r.ull);


    return 0;
}