SUBDIRS = burrus cwplib dfftpack djbfft-0.76 dsp dxml essl fftpack	\
fftw2 fxt goedecker gpfa green-ffts-2.0 gsl harm intel-mkl krukar	\
mixfft mpfun77 mpfun90 ooura qft sgimath singleton sorensen sunperf	\
temperton vdsp fftw3 athfft

EXTRA_DIST = Makefile.common

distclean-local:
	rm -f *.info *.data *.sysinfo

# this is a hack
RECURSIVE_TARGETS =
RECURSIVE_TARGETS += benchmark-recursive
benchmark: benchmark-recursive
benchmark-am:


