//BCH view, systematic encoding Reed Solomon using 3 bit symbols
#include "rs_gf8.h"
#include "gf8.h"

//#define OCTAL_MASK_A 		00707070707
//#define OCTAL_MASK_B 		07070707070
//#define BIT_PER_BYTE_MASK	0x0101010101010101


const uint32_t rs8_G_polys[] = {
	      01,	//0 symbols (dummy for indexing)
	     011,	//1 symbol
	    0132,	//2 symbols
	   01753,	//3 symbols
	  014775,	//4 symbols
	 0122313,	//5 symbols
	01576342	//6 symbols (identical to entries 1 thru 7 of the exponent table)
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
	for(uint32_t i = 1; i < raw; i <<= GF8_IDX_INC)
		msg_sz += GF8_IDX_INC;
	
	uint32_t chk = gf8_poly_mod(raw, msg_sz, chk_poly, chk_sz);
	raw <<= (chk_sz - GF8_IDX_INC);
	return raw | chk;
}

uint32_t rs8_get_syndromes(uint32_t p, int8_t p_sz, int8_t nsyms)
{
	uint32_t s = 0;
	for(--nsyms; nsyms >= 0; --nsyms)
	{
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
	for(int8_t i = 7; i > 0; --i)	//decrementing from 7 results in the inverse index so we don't need to call gf8_inverse(i)
	{
		if(e_pos & 1)
		{
			//faster equivalent of gf8_poly_mul() for a monic binomial
			e_loc = (e_loc << GF8_IDX_INC) ^ gf8_poly_scale(e_loc, gf8_exp[i]);
		}
		e_pos >>= 1;
	}

	return e_loc;
}

uint32_t rs8_get_errata_evaluator(uint32_t synd, uint32_t e_loc, int8_t chk_sz)
{
	uint32_t e_eval = gf8_poly_mul(synd, e_loc);
	return e_eval & ~((uint32_t)-1 << (chk_sz + GF8_IDX_INC));
}

uint32_t rs8_get_modified_syndromes(uint32_t synd, int8_t s_sz, uint32_t e_loc, int8_t e_sz)
{
	e_sz -= GF8_IDX_INC;
	//FIXME: under extreme conditions (check symbols + e_loc size -1 > 10), this can overflow
	// this can only happen at 6 check symbols so either a specialized version of poly_mul needs to be
	// made for here or we need strict limits on the number of errasures and errors that we attempt to correct
	// or simply don't support 6 check symbols since that would only leave 8 codewords

	// might be fine as it is though since it seems the terms we want are in the low end anyway

	uint32_t mod_synd = gf8_poly_mul(synd, e_loc) >> e_sz;
	mod_synd &= ~((uint32_t)-1 << (s_sz - e_sz));
	return mod_synd | (synd & ((uint32_t)-1 << e_sz));
}

/*
//Forney algorithm
uint32_t rs8_correct_errata(uint32_t recv, uint32_t synd, int8_t e_pos)
{

}*/