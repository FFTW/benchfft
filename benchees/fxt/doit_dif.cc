#define NAME "fxt-dif"
#define NOTES "Radix-4 DIF algorithm, complex datatype.  Backwards transform has explicit conjugation pass."
#define DOIT_FFT fft_dif4(x, m, is)

#include "doit.cc"
