#define NAME "arndt_fht"
#define NOTES "FFT by FHT algorithm, real/imag arrays datatype."

#define NO_Complex 1

#define DOIT_FFT fht_fft(inr, ini, m, is)

#include "doit.cc"
