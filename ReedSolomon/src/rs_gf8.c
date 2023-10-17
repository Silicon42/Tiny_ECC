//BCH view, systematic encoding Reed Solomon using 3 bit symbols
#include "rs_gf8.h"
#include "gf8.h"

//#define OCTAL_MASK_A 		00707070707
//#define OCTAL_MASK_B 		07070707070
//#define BIT_PER_BYTE_MASK	0x0101010101010101


const uint32_t rs8_G_polys[] = {
	      01,	//0 symbols (dummy for indexing)
	     012,	//1 symbol	First Consectutive Root, aka fcr aka c = 1
	    0163,	//2 symbols
	   01525,	//3 symbols
	  013123,	//4 symbols
	 0143562,	//5 symbols
	01111111	//6 symbols
};

/*
//returns a generator polynomial for BCH view Reed Solomon with a given number of check symbols
uint32_t rs8_G_unpack(int8_t syms)
{
	// packed octal format of the generator polynomials for given numbers of check symbols,
	// omits the trailing highest order 1 term 
	// indexing is done by shifting 3*n*(n-1)/2 where n is the number of check symbols
	// and then shifting and masking 3 bits at a time for a count of n
	uint64_t packed = 0576342223134775753321;
	
	polynomial coefficients for 1 thru 6 symbols highest to lowest term order are as follows
	          1,1	1 symbol
	        1,3,2	2 symbols
	      1,7,5,3	3 symbols
	    1,4,7,7,5	4 symbols
	  1,2,2,3,1,3	5 symbols
	1,5,7,6,3,4,2	6 symbols (identical to entries 1 thru 7 of the exponent table)
	
	int8_t chk_bits = syms * GF8_IDX_INC;
	packed >>= chk_bits * (syms-1) / 2;
	packed &= ~((uint32_t)-1 << chk_bits);
	return packed | (1 << chk_bits);
}
*/
//encodes a block of up to 18 bits worth of raw data as a Reed Solomon code word
//infers message length from provided data, doesn't verify that data length fits with the intended number of check symbols
//TODO: might be better to have the check polynomials easily indexable by check symbol count, consider rewriting for that
uint32_t rs8_encode(uint32_t raw, uint32_t chk_poly, int8_t chk_sz)
{
	int8_t msg_sz = 0;
	for(uint32_t i = 1; i <= raw; i <<= GF8_IDX_INC)
		msg_sz += GF8_IDX_INC;
	
	uint32_t chk = gf8_poly_mod(raw, msg_sz, chk_poly, chk_sz);
	raw <<= (chk_sz - GF8_IDX_INC);
	return raw | chk;
}

uint32_t rs8_get_syndromes(uint32_t p, int8_t p_sz, int8_t nsyms)
{
	uint32_t s = 0;	//accumulate syndromes in s
	for(; nsyms > 0; --nsyms)	//accumulation done in descending index order
	{	//which syndromes are used is effected by fcr so if you change that it must be changed here too
		s <<= GF8_IDX_INC;
		s |= gf8_poly_eval(p, p_sz, gf8_exp[nsyms]);
	}

	return s;
}

//e_pos is encoded such that a set bit indicates the corresponding degree term is erased or in error
// might not be faster than listing indices but is a bit more transparent and since the field is small
// should have relatively little impact.
//TODO: verify this ^
uint32_t rs8_get_errata_locator(int8_t e_pos)
{
	uint32_t e_loc = 1;	//the errata locator polynomial
	for(int8_t i = 0; i < 7; ++i)
	{
		if(e_pos & 1)
		{
			//faster equivalent of gf8_poly_mul() for a monic binomial in the form (ax - 1)
			//using this form of the roots simplifies later calculations since we know term 0 is a 1
			e_loc ^= gf8_poly_scale(e_loc, gf8_exp[i]) << GF8_IDX_INC;
		}
		e_pos >>= 1;
	}

	return e_loc;
}

uint32_t rs8_get_errata_evaluator(uint32_t synd, int8_t chk_sz, uint32_t e_loc)
{
	uint32_t e_eval;
	//TODO: this can be optimized since term 0 of e_loc is always 1 and no more than 5 terms beyond that are needed
	e_eval = gf8_poly_mul(synd, e_loc);
	return e_eval & ~((uint32_t)-1 << chk_sz);	//mask to the appropriate size
}
/*
uint32_t rs8_get_forney_syndromes(uint32_t synd, int8_t chk_sz, uint32_t e_loc, int8_t e_sz)
{
	e_sz -= GF8_IDX_INC;

	uint32_t fsynd = gf8_poly_mul(synd, e_loc);// >> e_sz;	//shift down by # of erasures worth of symbols
	uint32_t mask = (uint32_t)-1 << (chk_sz);// - e_sz);	//mask to keep the (# of check symbols - erasures) worth of symbols
	fsynd &= ~mask;
	return fsynd;// | (synd & mask);
}
*/
//Forney algorithm
uint32_t rs8_get_errata_magnitude(uint32_t synd, int8_t chk_sz, uint32_t e_loc, int8_t e_pos)
{
/*
	error value e(i) = -(X(i)^(1-c) * omega(X(i)^-1)) / (lambda'(X(i)^-1))
	where X(i)^-1 is the roots of the error locator, omega(X) is the error evaluator,
	lambda'(X) is the formal derivative of the error locator, and c is the 1st
	consecutive root of the generator used in the encoding, which in this case is always
	1 such that the X(i)^(1-c) simplifies out to 1

	e_eval = synd * e_loc
*/
	uint32_t e_eval, e_loc_prime, e_mag;
	int8_t root, ee_res, lp_res;

	e_eval = rs8_get_errata_evaluator(synd, chk_sz, e_loc);

	e_loc_prime = gf8_poly_formal_derivative(e_loc);

	e_mag = 0;
	for(int8_t i = 1; i < 8; ++i)
	{
		e_mag <<= GF8_IDX_INC;
		e_pos <<= 1;

		if(e_pos & 128)
		{
			root = gf8_exp[i];
			ee_res = gf8_poly_eval(e_eval, chk_sz, root);
			lp_res = gf8_poly_eval(e_loc_prime, chk_sz, root);	//chk_sz is guaranteed to be at least as big as e_loc_prime's actual size
			e_mag |= gf8_div(ee_res, lp_res);
		}
	}

	return e_mag;
}
/*
uint32_t rs8_correct_errata(uint32_t recv, uint32_t synd, int8_t e_pos)
{

}*/