#define NAME "fxt-dit"
#define NOTES "Radix-4 DIT algorithm, interleaved complex datatype.  Backwards transform has explicit conjugation pass."
#define DOIT_FFT fft_dit4(x, m, is)

#include "doit.cc"
