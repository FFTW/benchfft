if HAVE_F77
if DOUBLE
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

PKG = ffte-2.0
TARFILE = $(PKG).tgz
EXTRA_DIST = $(TARFILE) ffte-1.0.tgz zfft1d.f.patch zfft2d.f.patch zfft3d.f.patch
LOCAL_COPIES = cparam.h param.h

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@
doit.c: cparam.h

cparam.h: $(PKG)/param.h extract
	rm -f $@
	sed 's,^C\(.*\)$$,/* \1 */,;s,PARAMETER *(\([^)]*\)),static const int \1\;,' $(PKG)/param.h > $@

$(PKG)/param.h: extract

libfft_a_SOURCES = $(PKG)/kernel.f zfft1d.f zfft2d.f zfft3d.f
BUILT_SOURCES = zfft1d.f zfft2d.f zfft3d.f

$(PKG)/kernel.f $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f: extract

kernel.o zfft1d.o zfft2d.o zfft3d.o: param.h

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
