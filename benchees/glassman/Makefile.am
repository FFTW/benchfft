if HAVE_F77
if SINGLE
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

EXTRA_DIST = glassman-fft.f glassman-fft.patch doc
LOCAL_COPIES = glassman-fft-patched.f

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@

libfft_a_SOURCES = glassman-fft-patched.f
BUILT_SOURCES = glassman-fft-patched.f

glassman-fft-patched.f: glassman-fft.f glassman-fft.patch
	rm -f $@
	$(CP) $(srcdir)/glassman-fft.f $@
	$(PATCH) $@ <$(srcdir)/glassman-fft.patch

include ../Makefile.common
