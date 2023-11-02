// meant to benchmark some low level operations that get used frequently to determine which implementation to go with
#include "rs_gf8.h"
#include "gf8.h"
#include <stdio.h>
#include <stdint.h>
#include <time.h>

#define MILLION 1000000
#define BILLION 1000000000

int main(void)
{
	struct timespec start, end;
	int x = 1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	for (int i = 0; i < BILLION; i++)
	{
		// x ^= gf8_poly_scale(i, i);
		x ^= gf8_mul(i & 7, x & 7);
		x ^= gf8_mul(i & 7, x & 7);
		x ^= gf8_mul(i & 7, x & 7);
		x ^= gf8_mul(i & 7, x & 7);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
	printf("%o\n", x); // prevent optimizing x out;
	int64_t nanoseconds = (end.tv_sec - start.tv_sec) * BILLION;
	nanoseconds += (end.tv_nsec - start.tv_nsec);
	double seconds = (double)nanoseconds / (double)BILLION;
	printf("1 billion calls took: %f seconds\n", seconds);
	return 0;
}