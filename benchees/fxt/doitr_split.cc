#define NAME "fxt-split-real"
#define NOTES "FFT by real-data split-radix algorithm.  Original Fortran code by Sorensen; published in H.V. Sorensen et al., Real-valued fast Fourier transform algorithms, IEEE Trans. Acoust. Speech Sig. Proc. 35, 849-863 (1987).  Adapted to C by Bill Simpson, 1995.  Further optimizations by Joerg Arndt."

#define DOIT_FFT split_radix_real_complex_fft(in, m, -1)
#define DOIT_IFFT split_radix_complex_real_fft(in, m, +1)

#define HALFCOMPLEX 1

#include "doitr.cc"
