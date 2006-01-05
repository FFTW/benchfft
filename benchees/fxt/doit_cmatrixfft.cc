#define NAME "fxt-cmatrixfft"
#define NOTES "Four-step algorithm, uses FHT for sub-FFTs, interleaved real/imag format."

#define DOIT_FFT matrix_fft(x, m, is)

#include "doit.cc"
