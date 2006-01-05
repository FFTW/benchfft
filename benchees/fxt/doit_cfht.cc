#define NAME "fxt-cfht"
#define NOTES "FFT by FHT algorithm, interleaved real/imag format."

#define DOIT_FFT fht_fft(x, m, is)

#include "doit.cc"
