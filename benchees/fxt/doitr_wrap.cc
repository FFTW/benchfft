#define NAME "fxt-wrap-real"
#define NOTES "Real-data FFT by wrapping half-size complex FFT."

#define DOIT_FFT wrap_real_complex_fft(in, m, -1)
#define DOIT_IFFT wrap_complex_real_fft(in, m, +1)

#define PACKED 1

#include "doitr.cc"
