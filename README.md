# benchFFT

benchFFT is a program to benchmark FFT software, assembled by Matteo
Frigo and Steven G. Johnson at the Massachusetts Institute of
Technology as part of the FFTW project:

                http://www.fftw.org/benchfft

The benchmark incorporates a large number of publicly available FFT
implementations, in both C and Fortran, and measures their performance
and accuracy over a range of transform sizes. It benchmarks both real
and complex transforms in one, two, and three dimensions. The FFT implementations in the benchmark (found in the benchees/
subdirectory) were written by various authors over a period of more
than 35 years. Except for small tweaks to get things to compile, we
used the unmodified original codes.

# Warning

This code is rather old and may require some effort to compile on a modern machine.   (Run `bootstrap.sh` first on a git clone to auto-generate the build files.)
