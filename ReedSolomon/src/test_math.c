#include "rs_gf8.h"
#include <stdio.h>

int main()
{
    union Poly8 p, q, r;
    p.ull = 0x010307;
    q.ull = 0x305;

    r = gf8_poly_mul(p, q);

    printf("%llX", r.ull);

   // r = gf8_rs_G_unpack(3);

    //printf("\n%llX", r.ull);


    return 0;
}