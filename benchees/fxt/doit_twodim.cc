#define NAME "fxt-twodim"
#define NOTES "Uses FHT for 1d FFTs."

#define NO_Complex 1

#define RANK_OK(r) ((r) == 2)

#define DOIT_FFT twodim_fft(inr, ini, n0, n1, is)

#include "doit.cc"
