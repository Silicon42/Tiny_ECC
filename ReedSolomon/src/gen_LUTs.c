
#include <stdio.h>
#include "gf8.h"
#include "gf16.h"

int main ()
{
	int8_t x;
	int8_t exp_LUT[15];
	int8_t log_LUT[16];
	log_LUT[0] = 0xFF;
	uint64_t g_poly;
	
	printf("\nGenerating GF(16) LUTs\n");

	printf("exp LUT:\n");
	x = 1;
	for(int i = 0; i < 15; ++i)
	{
		printf("0x%X,", x);
		exp_LUT[i] = x;
		log_LUT[x] = i;
		x = gf16_mul2_noLUT(x);
	}
	printf("\nduplicate exp entries after to avoid modulo op\n");

	printf("log LUT:\n");
	for(int i = 0; i < 16; ++i)
	{
		printf("0x%X,", (uint8_t)log_LUT[i]);
	}

	printf("\nReed Solomon generator polynomials:\n0x1,\n");
	g_poly = 1;
	for(int i = 1; i < 15; ++i)
	{
		g_poly = gf16_poly_scale(g_poly, exp_LUT[i]) ^ (g_poly << GF16_SYM_SZ);
		printf("0x%llX,\n", g_poly);
	}


	printf("\n\nGenerating GF(8) LUTs\n");
	printf("exp LUT:\n");
	x = 1;
	for(int i = 0; i < 7; ++i)
	{
		printf("%i, ", x);
		exp_LUT[i] = x;
		log_LUT[x] = i;
		x = gf8_mul2_noLUT(x);
	}
	printf("\nduplicate exp entries after to avoid modulo op\n");

	printf("log LUT:\n");
	for(int i = 0; i < 8; ++i)
	{
		printf("%i, ", log_LUT[i]);
	}

	printf("\nReed Solomon generator polynomials:\n01,\n");
	g_poly = 1;
	for(int i = 1; i < 7; ++i)
	{
		g_poly = gf8_poly_scale(g_poly, exp_LUT[i]) ^ (g_poly << GF8_SYM_SZ);
		printf("0%o,\n", (uint32_t)g_poly);
	}

}