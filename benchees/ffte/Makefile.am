if HAVE_F77
if DOUBLE
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

PKG = ffte-1.0
TARFILE = $(PKG).tgz
EXTRA_DIST = $(TARFILE)
LOCAL_COPIES = param.h

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@
doit.c: param.h

param.h: $(PKG)/param.h extract
	rm -f $@
	sed 's,^C\(.*\)$$,/* \1 */,;s,PARAMETER *(\([^)]*\)),static const int \1\;,' $(PKG)/param.h > $@

$(PKG)/param.h: extract

libfft_a_SOURCES = $(PKG)/kernel.f $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f
$(PKG)/kernel.f $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f: extract

kernel.o zfft1d.o zfft2d.o zfft3d.o: $(PKG)/param.h

include ../Makefile.common
