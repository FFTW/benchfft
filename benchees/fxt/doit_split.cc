#define NAME "fxt-split"
#define NOTES "Split-radix algorithm, iterative trig. factor calculation."

#define NO_Complex 1

#define DOIT_FFT split_radix_fft(inr, ini, m, is)

#include "doit.cc"
