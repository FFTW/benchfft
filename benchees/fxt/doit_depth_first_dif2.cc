#define NAME "fxt-depth-first-dif2"
#define NOTES "Radix-2 DIF with swapped loop order, interleaved real/imag format."

#define DOIT_FFT fft_depth_first_dif2(x, m, is)

#include "doit.cc"
