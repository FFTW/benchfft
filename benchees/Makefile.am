SUBDIRS = athfft bloodworth burrus cwplib dfftpack djbfft-0.76 dsp dxml	\
emayer essl fftpack fftw2 fftw3 fxt goedecker gpfa green-ffts-2.0 gsl	\
harm intel-mkl krukar mfft mixfft mpfun77 mpfun90 nag napack nr ooura	\
qft rmayer sgimath singleton sorensen statlib sunperf temperton vdsp

EXTRA_DIST = Makefile.common

distclean-local:
	rm -f *.info *.data *.sysinfo

# this is a hack
RECURSIVE_TARGETS =
RECURSIVE_TARGETS += benchmark-recursive
benchmark: benchmark-recursive
benchmark-am:


