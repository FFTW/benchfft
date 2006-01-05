#define NAME "fxt-csplit"
#define NOTES "Split-radix algorithm, interleaved real/imag format."

#define DOIT_FFT split_radix_fft(x, m, is)

#include "doit.cc"
