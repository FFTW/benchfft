#define NAME "fxt-fht-real"
#define NOTES "FFT by FHT algorithm."

#define DOIT_FFT fht_real_complex_fft(in, m, -1)
#define DOIT_IFFT fht_complex_real_fft(in, m, +1)

#define HALFCOMPLEX 1

#include "doitr.cc"
