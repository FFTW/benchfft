SUBDIRS = acml arprec bloodworth burrus cross cwplib dfftpack dsp dxml	\
emayer esrfft essl ffmpeg ffte fftj fftpack fftreal fftw2 fftw3 fxt	\
glassman goedecker gpfa green-ffts-2.0 gsl harm hp-mlib imsl intel-mkl	\
intel-ipps jmfft kissfft krukar mfft mixfft monnier morris mpfun77	\
mpfun90 nag napack newsplit nr numutils ooura qft ransom rmayer		\
scimark2c sciport sgimath singleton sorensen spiral-fft statlib		\
sunperf temperton teneyck valkenburg vbigdsp vdsp

EXTRA_DIST = Makefile.common

distclean-local:
	rm -f *.info *.speed *.accuracy *.sysinfo

# handle recursive ``make benchmark'' etc.
BENCHFFT_RECURSIVE_TARGETS = benchmark-recursive accuracy-recursive	\
fftinfo-recursive

benchmark: benchmark-recursive
benchmark-am:
accuracy: accuracy-recursive
accuracy-am:
fftinfo: fftinfo-recursive
fftinfo-am:

$(BENCHFFT_RECURSIVE_TARGETS):
	@failcom='exit 1'; \
	for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=* | --[!k]*);; \
	    *k*) failcom='fail=yes';; \
	  esac; \
	done; \
	dot_seen=no; \
	target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    dot_seen=yes; \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	  || eval $$failcom; \
	done; \
	if test "$$dot_seen" = "no"; then \
	  $(MAKE) $(AM_MAKEFLAGS) "$$target-am" || exit 1; \
	fi; test -z "$$fail"
