SUBDIRS = burrus cwplib dfftpack djbfft-0.76 dsp dxml fftpack fftw2	\
fxt goedecker gpfa green-ffts-2.0 gsl harm intel-mkl krukar mixfft	\
mpfun77 mpfun90 ooura sgimath singleton sorensen sunperf temperton	\
vdsp

EXTRA_DIST = Makefile.common

distclean-local:
	rm -f *.info *.data *.sysinfo

# this is a hack
RECURSIVE_TARGETS =
RECURSIVE_TARGETS += benchmark-recursive
benchmark: benchmark-recursive
benchmark-am:


