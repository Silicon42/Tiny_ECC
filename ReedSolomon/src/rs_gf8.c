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
	uint32_t synd = 0;	//accumulate syndromes in s
	for(; nsyms > 0; --nsyms)	//accumulation done in descending index order
	{	//which syndromes are used is effected by fcr so if you change that it must be changed here too
		synd <<= GF8_SYM_SZ;
		synd |= gf8_poly_eval(p, p_sz, gf8_exp[nsyms]);
	}

	return synd;
}

//erase_pos is encoded such that a set bit indicates the corresponding degree term is erased or in error
// might not be faster than listing indices but is a bit more transparent and since the field is small
// should have relatively little impact.
//TODO: verify this ^
uint32_t rs8_get_erasure_locator(int8_t erase_pos)
{
	uint32_t erase_loc = 1;	//the errata locator polynomial
	for(int8_t i = 0; i < 7; ++i)
	{
		if(erase_pos & 1)
		{
			//faster equivalent of gf8_poly_mul() for a monic binomial in the form (ax - 1)
			//using this form of the roots simplifies later calculations since we know term 0 is a 1
			erase_loc ^= gf8_poly_scale(erase_loc, gf8_exp[i]) << GF8_SYM_SZ;
		}
		erase_pos >>= 1;
	}

	return erase_loc;
}

//this can also be used to get the Forney Syndromes
uint32_t rs8_get_errata_evaluator(uint32_t synd, gf8_idx chk_sz, uint32_t errata_loc)
{
	uint32_t errata_eval;
	//TODO: this can be optimized since term 0 of errata_loc is always 1 and no more than 5 terms beyond that are needed
	errata_eval = gf8_poly_mul(synd, errata_loc);
	errata_eval &= ~((uint32_t)-1 << chk_sz);	//mask to the appropriate size
	return errata_eval;
}

//Forney algorithm
uint32_t rs8_get_errata_magnitude(uint32_t errata_eval, gf8_idx chk_sz, uint32_t errata_loc, int8_t errata_pos)
{
/*
	error value e(i) = -(X(i)^(1-c) * omega(X(i)^-1)) / (lambda'(X(i)^-1))
	where X(i)^-1 is the roots of the error locator, omega(X) is the error evaluator,
	lambda'(X) is the formal derivative of the error locator, and c is the 1st
	consecutive root of the generator used in the encoding, which in this case is always
	1 such that the X(i)^(1-c) simplifies out to 1

	errata_eval = synd * errata_loc
*/
	uint32_t errata_loc_prime, errata_mag;
	int8_t root, ee_res, lp_res;

	errata_loc_prime = gf8_poly_formal_derivative(errata_loc);

	errata_mag = 0;
	for(int8_t i = 1; i < 8; ++i)
	{
		errata_mag <<= GF8_SYM_SZ;
		errata_pos <<= 1;

		if(errata_pos & 128)
		{
			root = gf8_exp[i];
			ee_res = gf8_poly_eval(errata_eval, chk_sz, root);
			lp_res = gf8_poly_eval(errata_loc_prime, chk_sz, root);	//chk_sz is guaranteed to be at least as big as errata_loc_prime's actual size
			errata_mag |= gf8_div(ee_res, lp_res);	//TODO: if converting to position list form, consider adding a pairwise divide function to gf8.c
		}
	}

	return errata_mag;
}

//synd_rem is the number of remaining syndromes, ie # check symbols - # erasures, aka N on Wikipedia
uint32_t rs8_get_error_locator(uint32_t synd, gf8_idx s_sz)
{
	uint32_t error_loc, error_loc_last, error_loc_temp;
	int8_t disc, disc_last;
	gf8_idx delay, error_sz;

	error_loc = 1;		//aka C(x)
	error_loc_last = 1;	//aka B(x)
	//error_loc_temp;	//aka T(x)
	//s_sz;				//aka N		(* GF8_SYM_SZ so it's in bits)
	error_sz = 0;		//aka L		(* GF8_SYM_SZ so it's in bits due to being involved in calculations with n)
	delay = GF8_SYM_SZ;	//aka m		(* GF8_SYM_SZ so it's in bits)
	disc_last = 1;		//aka b
	//disc;				//aka d

	for(gf8_idx n = 0; n < s_sz; n = gf8_idx_inc(n))
	{
		disc = (synd >> n) & 7;	//term 0 of the following pairwise product
		for(gf8_idx i = GF8_SYM_SZ; i <= error_sz; i = gf8_idx_inc(i))
		{	//TODO: consider adding a pairwise product function to gf8.c
			disc ^= gf8_mul((error_loc >> i) & 7, (synd >> (n - i)) & 7);
		}

		if(disc)
		{
			//this is only okay to have here because we're dealing with small fields, so we can make
			// use of automatic register renaming to not incur any significant extra cost from copying,
			// otherwise you would need the else block to avoid the copy as in the Wikipedia article
			error_loc_temp = error_loc;
			error_loc ^= (gf8_poly_scale(error_loc_last, gf8_div(disc, disc_last)) << delay);

			if(2 * error_sz <= n)
			{
				error_loc_last = error_loc_temp;
				error_sz = gf8_idx_inc(n) - error_sz;
				disc_last = disc;
				delay = 0;	//doesn't make sense to reset to 1 term of delay and continue since there's only 1 instruction before the end of the loop anyway
			}
		}
		delay = gf8_idx_inc(delay);
	}

	return error_loc;
}

//mask_pos has set bits for the valid (ie received) message terms and lets us skip the erasures and
// otherwise non-transmitted terms such as fixed padding or less than maximal message length
int8_t rs8_get_error_pos(uint32_t error_loc, int8_t mask_pos)
{
	int8_t error_pos = 0;

	for(int8_t i = 1; i <= 7; ++i)
	{
		error_pos <<= 1;
		mask_pos <<= 1;
		if(mask_pos & 0b10000000)	//skips non-received symbols, not strictly required but potentially beneficial since poly eval is relatively expensive
			error_pos |= !gf8_poly_eval(error_loc, 21, gf8_exp[i]);	//for non-C coders, this means that when it evaluates to 0 we get back a True which is equivalent to 1
	}

	return error_pos;
}

//intended r_sz of the received message must be specified b/c leading (ie high order) 0 terms are functionally the
// same as having a lower degree limit, however r_sz is used to determine valid error positions which can't occur
// in non-transmitted 0 padding, ie if 5 3-bit terms were transmitted
//tx_pos inludes set bits for only the valid positions for errors to occur, ie not in untransmitted padding symbols
uint32_t rs8_decode(uint32_t recv, gf8_idx r_sz, int8_t chk_syms, int8_t e_pos, int8_t tx_pos)
{
	int8_t erase_cnt = __builtin_popcount(e_pos);
	if(erase_cnt > chk_syms)	//if the number of erasures is greater than the number of check symbols,
		return -1;	// it's already beyond the Singleton Bound and can't be uniquely decoded so we return and flag an error

	uint32_t errata_eval = rs8_get_syndromes(recv, r_sz, chk_syms);

	if(errata_eval == 0)	//no errors
		return recv;

	gf8_idx chk_sz = chk_syms*GF8_SYM_SZ;
	uint32_t e_loc = 1;

	if(e_pos)	//this check isn't required but shortcuts excess calculations when no erasures specified
	{
		e_loc = rs8_get_erasure_locator(e_pos);
		errata_eval = rs8_get_errata_evaluator(errata_eval, chk_sz, e_loc);	//compute Forney syndromes
	}

	if(erase_cnt != chk_syms)	//skip checking for errors if the maximum number of erasures occurred as we no longer have enough extra data
	{	
		uint32_t error_loc = rs8_get_error_locator(errata_eval, chk_sz);	//may be smaller than chk_sz but under most conditions this is correct
		int8_t error_loc_order = gf8_poly_get_order(error_loc);
		if(2*error_loc_order > chk_syms - erase_cnt)	//check that the number of errors isn't beyond the Singleton Bound
			return -2;	//TODO: more diagnostic information should be encoded here and below
		int8_t error_pos = rs8_get_error_pos(error_loc, tx_pos & (~e_pos));
		int8_t error_cnt = __builtin_popcount(error_pos);
		if(error_cnt != error_loc_order)
			return -3;	//not enough or too many roots
		
		e_pos |= error_pos;
		e_loc = gf8_poly_mul(e_loc, error_loc);
		errata_eval = rs8_get_errata_evaluator(errata_eval, chk_sz, error_loc);
	}

	uint32_t errata_mag = rs8_get_errata_magnitude(errata_eval, chk_sz, e_loc, e_pos);
	recv ^= errata_mag;

	return recv;
}
