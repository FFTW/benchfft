if HAVE_F77
PRG=doit
if SINGLE
noinst_LIBRARIES=libffts.a
LIBFFT = libffts.a
else
noinst_LIBRARIES=libfftd.a
LIBFFT = libfftd.a
endif
endif

EXTRA_DIST = readme ffts.f fftd.f ffts.patch fftd.patch
LOCAL_COPIES = ffts-patched.f fftd-patched.f

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBFFT) $(LIBBENCH) @FLIBS@

libffts_a_SOURCES = ffts-patched.f
libfftd_a_SOURCES = fftd-patched.f
BUILT_SOURCES = ffts-patched.f fftd-patched.f

ffts-patched.f: ffts.f ffts.patch
	rm -f $@
	$(CP) $(srcdir)/ffts.f $@
	$(PATCH) $@ <$(srcdir)/ffts.patch

fftd-patched.f: fftd.f fftd.patch
	rm -f $@
	$(CP) $(srcdir)/fftd.f $@
	$(PATCH) $@ <$(srcdir)/fftd.patch

include ../Makefile.common
