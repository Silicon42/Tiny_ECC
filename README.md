# Tiny_ECC
Lightweight C implementations of short Reed-Solomon error correcting codes. Suitable for use on 32-bit embedded systems with no additional hardware or architecture support required to be performant.

Currently includes encoding for Reed-Solomon over GF(8) and GF(16) and decoding for Reed-Solomon over GF(8) with support for both erasures and errors.

GF(16) decoding coming soon.

May later add Binary Golay codes and Hamming codes if I ever have time for it.

## Reason for writing
I needed some very small block length error correcting codes suitable for use in designing fiducial markers and being used on embedded systems but couldn't find any existing open source ones that fit my needs so I decided I'd learn how it worked and do it myself.

I also realized that for small fields, roughly GF(2^7) or smaller depending on native register sizes, with a careful implementation, they could potentially compute certain operations in parallel at a symbol level without any additional architecture support and perform faster than the traditional method of iterating over terms using logarithm and exponent look up tables.

Somewhere along the line it also turned into an exercise in trying to write a more understandable toy example for anyone else trying to learn how to implement Reed-Solomon because for something that's used in as many places as it is today, the explanations typically are either overly simplified or are steeped in post-graduate level math concepts. They are also almost entirely focussed on GF(2^8), which while probably the most practically useful, is hell to write tests for because doing operations in GF(2^8) by hand to verify that a calculation is correct is extremely tedious and potentially error prone. From my perspective a better alternative is to do the same over GF(2^3) and then scale it up once all the conceptual kinks are ironed out.

## Sources I used for learning
Wikiversity:
+ [Reed–Solomon codes for coders](https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders)

I primarily used this and the included example as reference and was mostly satisfied with it but the unnecessarily and totally unexplained reversed order of the sense of some of the polynomials tripped me up. Parts of the explanation were rather lacking and there were a couple mistakes in it that were counteracted in other parts of it down the line but it was at least approachable as opposed to having loads of abstract math in matrix form thrown at me.

Wikipedia:
+ [Reed–Solomon error correction](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction)
+ [Forney algorithm](https://en.wikipedia.org/wiki/Forney_algorithm)
+ [Berlekamp–Massey algorithm](https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm)

Used as a second reference point to clear my confusion about the Wikiversity example. Very math notation heavy but makes a bit more sense once you can connect actual meaningful names to the variables as opposed to arbitrary Greek letters.

libcorrect:
+ [GitHub](https://github.com/quiet/libcorrect):

Third reference point used to understand how the erasure and error correction process mechanically works. Cleaner implementation in C++, not perfect because some of the function calls are a bit opaque since they just pass a reference to the entire Reed-Solomon data structure instead of the relevant bits of it.

_Algebraic Codes for Data Transmission_ by Richard E. Blahut
+ Chapter 7.4 Decoding of nonbinary BCH codes

When all else fails go back to the source. I still had some confusion about how it should work when the above examples seemed to contradict each other and this cleared it right up.
