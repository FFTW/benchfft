#define NAME "fxt-matrixfft"
#define NOTES "Four-step algorithm, real/imag arrays datatype, uses FHT for sub-FFTs."

#define NO_Complex 1

#define DOIT_FFT matrix_fft(inr, ini, m, is)

#include "doit.cc"
