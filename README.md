# Tiny_ECC
Lightweight C implementations of short error correcting codes


Currently includes encoding for Reed-Solomon over GF(8) and GF(16).

Decoding not yet implemented.

May later add Binary Golay encoding and Decoding and Hamming if I ever get around to it.

## Reason for writing
I needed some very small block length error correcting codes but couldn't find any existing open source ones that fit my needs so I decided I'd learn how it worked and do it myself.

## Sources I used for learning
[Wikiversity](https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders): I primarily used this and the included example as reference and was mostly satisfied with it but the unnecessarily and totally unexplained reversed order of the sense of some of the polynomials seriously tripped me up.
[Wikipedia](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction): Used as a second reference point to clear my confusion about the Wikiversity example.