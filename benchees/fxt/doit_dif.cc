#define NAME "arndt_dif"
#define NOTES "Radix-4 DIF algorithm, complex datatype.  Backwards transform has explicit conjugation pass."
#define DOIT_FFT dif4_fft(x, m, is)

#include "doit.cc"
