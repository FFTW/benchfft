SUBDIRS = arprec athfft bloodworth burrus cross cwplib dfftpack		\
djbfft-0.76 dsp dxml emayer esrfft essl ffte fftpack fftw2 fftw3 fxt	\
glassman goedecker gpfa green-ffts-2.0 gsl harm imsl intel-mkl jmfft	\
krukar mfft mixfft monnier mpfun77 mpfun90 nag napack nr numutils	\
ooura pfftw qft ransom rmayer sciport sgimath singleton sorensen	\
spiral-fft statlib sunperf temperton teneyck valkenburg vbigdsp vdsp

EXTRA_DIST = Makefile.common

distclean-local:
	rm -f *.info *.speed *.accuracy *.sysinfo

# this is a hack
RECURSIVE_TARGETS =
RECURSIVE_TARGETS += benchmark-recursive accuracy-recursive
benchmark: benchmark-recursive
benchmark-am:
accuracy: accuracy-recursive
accuracy-am:


