#define NAME "fxt-dif4"
#define NOTES "Radix-4 DIF algorithm, split real/imag format"

#define NO_Complex 1

#define DOIT_FFT fft_dif4(inr, ini, m, is)

#include "doit.cc"
