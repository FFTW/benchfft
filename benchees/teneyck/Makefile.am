if HAVE_F77
if SINGLE
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

EXTRA_DIST = fftlib.f fftlib.patch doc
LOCAL_COPIES = fftlib-patched.f

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@

libfft_a_SOURCES = fftlib-patched.f
BUILT_SOURCES = fftlib-patched.f

fftlib-patched.f: fftlib.f fftlib.patch
	rm -f $@
	$(CP) $(srcdir)/fftlib.f $@
	$(PATCH) $@ <$(srcdir)/fftlib.patch


include ../Makefile.common
