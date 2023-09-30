
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "rs_gf8.h"

int main ()
{
	uint8_t x = 1;

	/*
	printf("Generating GF(2^4) multiply LUT\n\n");

	for (int i = 0; i < 16; i++)
	{
		for(int j = 0; j < 16; j++)
		{
			printf("0x%X,", cl_mul4_noLUT(i, j));
		}
		printf("\n");
	}

	printf("\nGenerating GF(2^4) log & exp LUT\n\n");
	uint8_t gf2_4_log_LUT [15];
	uint8_t gf2_4_exp_LUT [16];

	for (int i = 0; i < 15; i++)
	{
		gf2_4_exp_LUT[i] = x;
		gf2_4_log_LUT[x] = i;
		x = cl_mul4_noLUT(x, 2);
	}

	printf("\nexp LUT\n");
	for (int i = 0; i < 15; i++)
	{
		printf("0x%X,", gf2_4_exp_LUT[i]);
	}
	printf("\nduplicate the values after to avoid modulo op\n");

	printf("\nlog LUT\n");
	for (int i = 0; i < 16; i++)
	{
		printf("0x%X,", gf2_4_log_LUT[i]);
	}


	printf("\n\nGenerating GF(2^3) multiply LUT\n\n");

	for (int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 8; j++)
		{
			printf("0x%X,", cl_mul3_noLUT(i, j));
		}
		printf("\n");
	}
*/
	printf("\nGenerating GF(2^3) log & exp LUT\n\n");
	uint8_t gf8_log_LUT [8] = {0xF8};
	uint8_t gf8_exp_LUT [7];

	for (int i = 0; i < 7; i++)
	{
		//printf("%i, %i\n",i,x);
		gf8_exp_LUT[i] = x;
		gf8_log_LUT[x] = i;
		x = gf8_mul2_noLUT(x);
	}

	printf("\nexp LUT\n");
	for (int i = 0; i < 7; i++)
	{
		printf("%X, ", gf8_exp_LUT[i]);
	}
	printf("\nduplicate the values after to avoid modulo op\n");

	printf("\nlog LUT\n0x");
	for (int i = 0; i < 8; i++)
	{
		printf("%X, ", gf8_log_LUT[i]);
	}

}