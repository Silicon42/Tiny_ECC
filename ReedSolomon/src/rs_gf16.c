// BCH view, systematic encoding Reed Solomon using 3 bit symbols
#include "rs_gf16.h"

#define RS16_BLOCK_MASK 0xFFFFFFFFFFFFFFF // mask that represents the valid symbol positions

const gf16_poly rs16_G_polys[] = {
				  0x1,	// 0 symbols (dummy for indexing)
				 0x12,	// 1 symbol	First Consectutive Root, aka fcr aka c = 1
				0x168,	// 2 symbols
			   0x1EDC,	// 3 symbols
			  0x1DC87,	// 4 symbols
			 0x1B4621,	// 5 symbols
			0x1793CAC,	// 6 symbols
		   0x1CDF27ED,	// 7 symbols
		  0x19434D6EC,	// 8 symbols
		 0x131D93D7A1,	// 9 symbols
		0x148AC942C27,	// 10 symbols
	   0x1A53AD3F368C,	// 11 symbols
	  0x1595814D94CD8,	// 12 symbols
	 0x185A439C7BDE62,	// 13 symbols
	0x111111111111111	// 14 symbols
};

// encodes a block of up to 56 bits worth of raw data as a Reed Solomon code word
//  infers message length from provided data, doesn't verify that data length fits
//  with the specified number of check symbols and will truncate to the max size,
//  discarding the most significant bits if oversized.
gf16_poly rs16_encode(gf16_poly raw, int8_t chk_syms)
{
	gf16_idx chk_sz = chk_syms * GF16_SYM_SZ;
	raw &= RS16_BLOCK_MASK >> chk_sz;	//truncate most significant bits if provided data is oversized
	gf16_idx msg_sz = gf16_poly_get_size(raw);
	chk_sz += GF16_SYM_SZ;

	gf16_poly chk = gf16_poly_mod(raw, msg_sz, rs16_G_polys[chk_syms], chk_sz);
	raw <<= chk_sz - GF16_SYM_SZ;
	return raw | chk;
}

gf16_poly rs16_get_syndromes(gf16_poly p, gf16_idx p_sz, int8_t nsyms)
{
	gf16_poly synd = 0;	// accumulate syndromes in synd
	for (; nsyms > 0; --nsyms)	// accumulation done in descending index order
	{	// which syndromes are used is effected by fcr so if you change that it must be changed here too
		synd <<= GF16_SYM_SZ;
		synd |= gf16_poly_eval(p, p_sz, gf16_exp[nsyms]);
	}

	return synd;
}

// erase_pos is encoded such that a set bit indicates the corresponding degree term is erased or in error
//  might not be faster than listing indices but is a bit more transparent and since the field is small
//  should have relatively little impact.
// TODO: verify this ^
gf16_poly rs16_get_erasure_locator(int16_t erase_pos)
{
	gf16_poly erase_loc = 1;	// the errata locator polynomial
	for (int8_t i = 0; i < GF16_MAX; ++i)
	{
		if (erase_pos & 1)
		{
			// faster equivalent of gf16_poly_mul() for a monic binomial in the form (ax - 1)
			// using this form of the roots simplifies later calculations since we know term 0 is a 1
			erase_loc ^= gf16_poly_scale(erase_loc, gf16_exp[i]) << GF16_SYM_SZ;
		}
		erase_pos >>= 1;
	}

	return erase_loc;
}

// this can also be used to get the Forney Syndromes
gf16_poly rs16_get_errata_evaluator(gf16_poly synd, gf16_idx chk_sz, gf16_poly errata_loc)
{
	// term 0 of errata_loc is always 1 and no more than 5 terms beyond that are needed
	gf16_poly errata_eval = gf16_poly_mul_q0_monic(synd, errata_loc);
	errata_eval &= ~((gf16_poly)-1 << chk_sz);	// mask to the appropriate size
	return errata_eval;
}

// Forney algorithm
gf16_poly rs16_get_errata_magnitude(gf16_poly errata_eval, gf16_idx chk_sz, gf16_poly errata_loc, int16_t errata_pos)
{
	/*
	error value e(i) = -(X(i)^(1-c) * omega(X(i)^-1)) / (lambda'(X(i)^-1))
	where X(i)^-1 is the roots of the error locator, omega(X) is the error evaluator,
	lambda'(X) is the formal derivative of the error locator, and c is the 1st
	consecutive root of the generator used in the encoding, which in this case is always
	1 such that the X(i)^(1-c) simplifies out to 1

	errata_eval = synd * errata_loc
	*/
	gf16_poly errata_loc_prime, errata_mag;
	gf16_elem root, ee_res, lp_res;

	errata_loc_prime = gf16_poly_formal_derivative(errata_loc);

	errata_mag = 0;
	for (int8_t i = 1; i <= GF16_MAX; ++i)
	{
		errata_mag <<= GF16_SYM_SZ;
		errata_pos <<= 1;

		if (errata_pos & 0x8000)
		{
			root = gf16_exp[i];
			ee_res = gf16_poly_eval(errata_eval, chk_sz, root);
			lp_res = gf16_poly_eval(errata_loc_prime, chk_sz, root);	// chk_sz is guaranteed to be at least as big as errata_loc_prime's actual size
			errata_mag |= gf16_div(ee_res, lp_res);	// TODO: if converting to position list form, consider adding a pairwise divide function to gf16.c
		}
	}

	return errata_mag;
}

// synd_rem is the number of remaining syndromes, ie # check symbols - # erasures, aka N on Wikipedia
gf16_poly rs16_get_error_locator(gf16_poly synd, gf16_idx s_sz)
{
	gf16_poly error_loc, error_loc_last, error_loc_temp;
	gf16_elem disc, disc_last;
	gf16_idx delay, error_sz;

	error_loc = 1;			// aka C(x)
	error_loc_last = 1;		// aka B(x)
	// error_loc_temp;		// aka T(x)
	// s_sz;				// aka N	(*GF16_SYM_SZ so it's in bits)
	error_sz = 0;			// aka L	(*GF16_SYM_SZ so it's in bits due to being involved in calculations with n)
	delay = GF16_SYM_SZ;	// aka m	(*GF16_SYM_SZ so it's in bits)
	disc_last = 1;			// aka b
	// disc;				// aka d

	for (gf16_idx n = 0; n < s_sz; n += GF16_SYM_SZ)
	{
		disc = (synd >> n) & GF16_MAX;	// term 0 of the following pairwise product
		for (gf16_idx i = GF16_SYM_SZ; i <= error_sz; i += GF16_SYM_SZ)
		{
			disc ^= gf16_mul((error_loc >> i) & GF16_MAX, (synd >> (n - i)) & GF16_MAX);
		}

		if (disc)
		{
			// this is only okay to have here because we're dealing with small fields, so we can make
			//  use of automatic register renaming to not incur any significant extra cost from copying,
			//  otherwise you would need the else block to avoid the copy as in the Wikipedia article
			error_loc_temp = error_loc;
			error_loc ^= (gf16_poly_scale(error_loc_last, gf16_div(disc, disc_last)) << delay);

			if (2 * error_sz <= n)
			{
				error_loc_last = error_loc_temp;
				error_sz = GF16_SYM_SZ + n - error_sz;
				disc_last = disc;
				delay = 0;	// doesn't make sense to reset to 1 term of delay and continue since there's only 1 instruction before the end of the loop anyway
			}
		}
		delay += GF16_SYM_SZ;
	}

	return error_loc;
}

// mask_pos has set bits for the valid (ie received) message terms and lets us skip the erasures and
//  otherwise non-transmitted terms such as fixed padding or less than maximal message length
int16_t rs16_get_error_pos(gf16_poly error_loc, int16_t mask_pos)
{
	int16_t error_pos = 0;

	for (int8_t i = 1; i <= GF16_MAX; ++i)
	{
		error_pos <<= 1;
		mask_pos <<= 1;
		if (mask_pos & 0b10000000)	// skips non-received symbols, not strictly required but potentially beneficial since poly eval is relatively expensive
			error_pos |= !gf16_poly_eval(error_loc, 21, gf16_exp[i]);	// for non-C coders, this means that when it evaluates to 0 we get back a True which is equivalent to 1
	}

	return error_pos;
}

// tx_pos inludes set bits for only the valid positions for errors to occur, ie not in untransmitted padding symbols
gf16_poly rs16_decode(gf16_poly recv, gf16_idx r_sz, int8_t chk_syms, int16_t e_pos, int16_t tx_pos)
{
	int8_t erase_cnt = __builtin_popcount(e_pos);
	if (erase_cnt > chk_syms)	// if the number of erasures is greater than the number of check symbols,
		return -1;	// it's already beyond the Singleton Bound and can't be uniquely decoded so we return an error value

	gf16_poly e_eval = rs16_get_syndromes(recv, r_sz, chk_syms);

	if (e_eval == 0) // no errors
		return recv;

	gf16_idx chk_sz = chk_syms * GF16_SYM_SZ;
	gf16_poly e_loc = 1;

	if (e_pos) // this check isn't required but shortcuts excess calculations when no erasures specified
	{
		e_loc = rs16_get_erasure_locator(e_pos);
		e_eval = rs16_get_errata_evaluator(e_eval, chk_sz, e_loc);	// compute Forney syndromes
	}

	if (erase_cnt != chk_syms)	// skip checking for errors if the maximum number of erasures occurred as we no longer have enough extra data
	{
		gf16_poly error_loc = rs16_get_error_locator(e_eval, chk_sz);	// may be smaller than chk_sz but under most conditions this is correct
		int8_t error_loc_order = gf16_poly_get_order(error_loc);
		if (2 * error_loc_order > chk_syms - erase_cnt)	// check that the number of errors isn't beyond the Singleton Bound
			return -2;	// TODO: more diagnostic information should be encoded here and below
		int16_t error_pos = rs16_get_error_pos(error_loc, tx_pos & (~e_pos));
		int8_t error_cnt = __builtin_popcount(error_pos);
		if (error_cnt != error_loc_order)
			return -3;	// not enough or too many roots

		// combine the error and erasure position, locator, and evaluator to the errata versions of themselves
		e_pos |= error_pos;
		e_loc = gf16_poly_mul(e_loc, error_loc);
		e_eval = rs16_get_errata_evaluator(e_eval, chk_sz, error_loc);
	}

	gf16_poly errata_mag = rs16_get_errata_magnitude(e_eval, chk_sz, e_loc, e_pos);
	recv ^= errata_mag;

	return recv;
}
