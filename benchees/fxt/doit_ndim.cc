#define NAME "fxt-ndim"
#define NOTES "Uses FHT for 1d FFTs."

#define NO_Complex 1

#define RANK_OK(r) ((r) >= 2 && (r) <= 5)

#define DOIT_FFT ndim_fft(inr, ini, rank, md, is)

#include "doit.cc"
