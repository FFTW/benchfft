#define NAME "fxt-cdif4"
#define NOTES "Radix-4 DIF algorithm, interleaved real/imag format"

#define DOIT_FFT fft_dif4(x, m, is)

#include "doit.cc"
