#define NAME "fxt-wrap-real"
#define NOTES "Real-data FFT by wrapping half-size complex FFT (via FHT)."

#define DOIT_FFT wrap_real_complex_fft(in, m)
#define DOIT_IFFT wrap_complex_real_fft(in, m)

#define PACKED 1

#include "doitr.cc"
