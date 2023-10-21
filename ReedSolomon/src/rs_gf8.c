//BCH view, systematic encoding Reed Solomon using 3 bit symbols
#include "rs_gf8.h"

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
	
	int8_t chk_bits = syms * GF8_SYM_SZ;
	packed >>= chk_bits * (syms-1) / 2;
	packed &= ~((uint32_t)-1 << chk_bits);
	return packed | (1 << chk_bits);
}
*/
//encodes a block of up to 18 bits worth of raw data as a Reed Solomon code word
//infers message length from provided data, doesn't verify that data length fits with the intended number of check symbols
//TODO: might be better to have the check polynomials easily indexable by check symbol count, consider rewriting for that
uint32_t rs8_encode(uint32_t raw, uint32_t chk_poly, gf8_idx chk_sz)
{
	gf8_idx msg_sz = 0;	//TODO: write a gf8_poly_get_size function to be used here
	for(uint32_t i = 1; i <= raw; i <<= GF8_SYM_SZ)
		msg_sz = gf8_idx_inc(msg_sz);
	
	uint32_t chk = gf8_poly_mod(raw, msg_sz, chk_poly, chk_sz);
	raw <<= gf8_idx_dec(chk_sz);
	return raw | chk;	//TODO: add minimal input protection to make sure output is within valid block length via a mask
}

uint32_t rs8_get_syndromes(uint32_t p, gf8_idx p_sz, int8_t nsyms)
{
	uint32_t s = 0;	//accumulate syndromes in s
	for(; nsyms > 0; --nsyms)	//accumulation done in descending index order
	{	//which syndromes are used is effected by fcr so if you change that it must be changed here too
		s <<= GF8_SYM_SZ;
		s |= gf8_poly_eval(p, p_sz, gf8_exp[nsyms]);
	}

	return s;
}

//e_pos is encoded such that a set bit indicates the corresponding degree term is erased or in error
// might not be faster than listing indices but is a bit more transparent and since the field is small
// should have relatively little impact.
//TODO: verify this ^
uint32_t rs8_get_erasure_locator(int8_t e_pos)
{
	uint32_t e_loc = 1;	//the errata locator polynomial
	for(int8_t i = 0; i < 7; ++i)
	{
		if(e_pos & 1)
		{
			//faster equivalent of gf8_poly_mul() for a monic binomial in the form (ax - 1)
			//using this form of the roots simplifies later calculations since we know term 0 is a 1
			e_loc ^= gf8_poly_scale(e_loc, gf8_exp[i]) << GF8_SYM_SZ;
		}
		e_pos >>= 1;
	}

	return e_loc;
}

//this can also be used to get the Forney Syndromes
uint32_t rs8_get_errata_evaluator(uint32_t synd, gf8_idx chk_sz, uint32_t e_loc)
{
	uint32_t e_eval;
	//TODO: this can be optimized since term 0 of e_loc is always 1 and no more than 5 terms beyond that are needed
	e_eval = gf8_poly_mul(synd, e_loc);
	return e_eval & ~((uint32_t)-1 << chk_sz);	//mask to the appropriate size
}

//Forney algorithm
uint32_t rs8_get_errata_magnitude(uint32_t e_eval, gf8_idx chk_sz, uint32_t e_loc, int8_t e_pos)
{
/*
	error value e(i) = -(X(i)^(1-c) * omega(X(i)^-1)) / (lambda'(X(i)^-1))
	where X(i)^-1 is the roots of the error locator, omega(X) is the error evaluator,
	lambda'(X) is the formal derivative of the error locator, and c is the 1st
	consecutive root of the generator used in the encoding, which in this case is always
	1 such that the X(i)^(1-c) simplifies out to 1

	e_eval = synd * e_loc
*/
	uint32_t e_loc_prime, e_mag;
	int8_t root, ee_res, lp_res;

	e_loc_prime = gf8_poly_formal_derivative(e_loc);

	e_mag = 0;
	for(int8_t i = 1; i < 8; ++i)
	{
		e_mag <<= GF8_SYM_SZ;
		e_pos <<= 1;

		if(e_pos & 128)
		{
			root = gf8_exp[i];
			ee_res = gf8_poly_eval(e_eval, chk_sz, root);
			lp_res = gf8_poly_eval(e_loc_prime, chk_sz, root);	//chk_sz is guaranteed to be at least as big as e_loc_prime's actual size
			e_mag |= gf8_div(ee_res, lp_res);	//TODO: if converting to position list form, consider adding a pairwise divide function to gf8.c
		}
	}

	return e_mag;
}

//synd_rem is the number of remaining syndromes, ie # check symbols - # erasures, aka N on Wikipedia
uint32_t rs8_get_error_locator(uint32_t synd, gf8_idx s_sz)
{
	uint32_t err_loc, err_loc_last, err_loc_temp;
	int8_t disc, disc_last;
	gf8_idx delay, err_cnt;

	err_loc = 1;		//aka C(x)
	err_loc_last = 1;	//aka B(x)
	//err_loc_temp;		//aka T(x)
	//s_sz;				//aka N		(* GF8_SYM_SZ so it's in bits)
	err_cnt = 0;		//aka L		(* GF8_SYM_SZ so it's in bits due to being involved in calculations with n)
	delay = GF8_SYM_SZ;	//aka m		(* GF8_SYM_SZ so it's in bits)
	disc_last = 1;		//aka b
	//disc;				//aka d

	for(gf8_idx n = 0; n < s_sz; n = gf8_idx_inc(n))
	{
		disc = (synd >> n) & 7;	//term 0 of the following pairwise product
		for(gf8_idx i = GF8_SYM_SZ; i < err_cnt; i = gf8_idx_inc(i))
		{	//TODO: consider adding a pairwise product function to gf8.c
			disc ^= gf8_mul((err_loc >> i) & 7, (synd >> (n - i)) & 7);
		}

		if(disc)
		{
			//this is only okay to have here because we're dealing with small fields,
			// so we can make use of automatic register renaming,
			// otherwise you would need the else block to avoid the copy as in the Wikipedia article
			err_loc_temp = err_loc;

			err_loc ^= (gf8_poly_scale(err_loc_last, gf8_div(disc, disc_last)) << delay);
			if(2 * err_cnt <= n)
			{
				err_loc_last = err_loc_temp;
				err_cnt = gf8_idx_inc(n) - err_cnt;
				disc_last = disc;
				delay = 0;	//doesn't make sense to reset to 1 term of delay and continue since there's only 1 instruction before the end of the loop anyway
			}
		}
		delay = gf8_idx_inc(delay);
	}

	return err_loc;
}

int8_t rs8_get_error_pos(uint32_t error_loc, int8_t erase_pos)
{
	int8_t error_pos = 0;
	int8_t bit = 1;
	erase_pos = ~erase_pos;	//simplifies the conditional experession

	for(int8_t i = 7; i > 0; --i)
	{
		if(erase_pos & bit)	//skips erasure positions, not strictly required but 
			error_pos |= gf8_poly_eval(error_loc, 21, i) ? 0 : bit;	//TODO: see if size calculation of error_loc would improve this over just taking the max size
			//^ marks the position of an error at a given bit index if the inverse of the index evaluates to 0
		bit <<= 1;
	}

	return error_pos;
}

//intended r_sz of the received message must be specified b/c leading (ie high order) 0 terms are functionally the
// same as having a lower degree limit, however r_sz is used to determine valid error positions which can't occur
// in non-transmitted 0 padding, ie if 5 3-bit terms were transmitted
uint32_t rs8_decode(uint32_t recv, gf8_idx r_sz, int8_t chk_syms, int8_t e_pos)
{
	int8_t erase_cnt = __builtin_popcount(e_pos);
	if(erase_cnt > chk_syms)	//if the number of erasures is greater than the number of check symbols,
		return -1;	// it's already beyond the Singleton Bound and can't be uniquely decoded so we return and flag an error

	uint32_t synd, e_loc, e_eval;

	synd = rs8_get_syndromes(recv, r_sz, chk_syms);

	if(synd == 0)	//no errors
		return recv;

	gf8_idx chk_sz = chk_syms*GF8_SYM_SZ;

	if(e_pos)
	{
		e_loc = rs8_get_erasure_locator(e_pos);
		synd = rs8_get_errata_evaluator(synd, chk_sz, e_loc);
	}

	if(synd != 0)	//errors remaining
	{
		uint32_t err_loc = rs8_get_error_locator(synd, chk_sz);	//FIXME: this should not be chk_sz but for testing it will do
		int8_t err_loc_order = gf8_poly_get_order(err_loc);
		if(2*err_loc_order > chk_syms - erase_cnt)	//check that the number of errors isn't beyond the Singleton Bound
			return -2;	//TODO: more diagnostic information should be encoded here and below
		int8_t err_pos = rs8_get_error_pos(err_loc, e_pos);
		int8_t err_cnt = __builtin_popcount(err_pos);
		if(err_cnt != err_loc_order)
			return -3;	//not enough or too many roots
		e_loc = gf8_poly_mul(e_loc, err_loc);
		e_eval = rs8_get_errata_evaluator(synd, chk_sz, err_loc);
		
	}
	
	return recv ^ rs8_get_errata_magnitude(e_eval, chk_sz, e_loc, e_pos);
}
