#define NAME "fxt-dit4"
#define NOTES "Radix-4 DIT algorithm, split real/imag format"

#define NO_Complex 1

#define DOIT_FFT fft_dit4(inr, ini, m, is)

#include "doit.cc"
