#define NAME "fxt-split"
#define NOTES "Split-radix algorithm, split real/imag format."

#define NO_Complex 1

#define DOIT_FFT split_radix_fft(inr, ini, m, is)

#include "doit.cc"
