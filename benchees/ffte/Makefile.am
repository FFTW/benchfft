if HAVE_F77
if DOUBLE
if HAVE_SSE2
PRG=doit doit_sse2
noinst_LIBRARIES=libfft.a libfft_kernel2.a libfft_sse2.a
else
PRG=doit
noinst_LIBRARIES=libfft.a libfft_kernel2.a
endif
endif
endif

PKG = ffte-3.0
TARFILE = $(PKG).tgz
EXTRA_DIST = $(TARFILE) zfft1d.f.patch zfft2d.f.patch zfft3d.f.patch
LOCAL_COPIES = cparam.h param.h

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a libfft_kernel2.a $(LIBBENCH) @FLIBS@
doit.c: doit-ffte.c cparam.h

doit_sse2_SOURCES=doit_sse2.c
doit_sse2_LDADD=libfft.a libfft_sse2.a $(LIBBENCH) @FLIBS@
doit_sse2.c: doit-ffte.c cparam.h

cparam.h: $(PKG)/param.h extract
	rm -f $@
	sed 's,^C\(.*\)$$,/* \1 */,;s,PARAMETER *(\([^)]*\)),static const int \1\;,' $(PKG)/param.h > $@

$(PKG)/param.h: extract

libfft_a_SOURCES = $(PKG)/kernel1.f zfft1d.f zfft2d.f zfft3d.f
BUILT_SOURCES = zfft1d.f zfft2d.f zfft3d.f

libfft_kernel2_a_SOURCES = $(PKG)/kernel2.f

libfft_sse2_a_SOURCES = $(PKG)/sse2.c
libfft_sse2_a_CFLAGS = $(AM_CFLAGS) @SSE2_CFLAGS@

$(PKG)/kernel1.f $(PKG)/kernel2.f $(PKG)/sse2.c $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f: extract

zfft1d.o zfft2d.o zfft3d.o: param.h

# Patch source code so that it compiles under F90 (there is a missing
# type declaration in one place).

zfft1d.f: $(PKG)/zfft1d.f zfft1d.f.patch
	rm -f $@
	$(CP) $(PKG)/zfft1d.f $@
	$(PATCH) $@ <$(srcdir)/zfft1d.f.patch

zfft2d.f: $(PKG)/zfft2d.f zfft2d.f.patch
	rm -f $@
	$(CP) $(PKG)/zfft2d.f $@
	$(PATCH) $@ <$(srcdir)/zfft2d.f.patch

zfft3d.f: $(PKG)/zfft3d.f zfft3d.f.patch
	rm -f $@
	$(CP) $(PKG)/zfft3d.f $@
	$(PATCH) $@ <$(srcdir)/zfft3d.f.patch

# need a local copy or Fortran include statement is confused
param.h: $(PKG)/param.h
	rm -f $@
	$(CP) $(PKG)/param.h $@

include ../Makefile.common
