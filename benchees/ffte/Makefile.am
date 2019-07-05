if HAVE_F77
if DOUBLE
PRG=doit doitv
noinst_LIBRARIES=libfft.a libfftv.a libfft_kernel.a
endif
endif

PKG = ffte-6.0
TARFILE = $(PKG).tgz
EXTRA_DIST = $(TARFILE) vzfft2d.f.patch vzfft3d.f.patch doit-ffte.c
LOCAL_COPIES = cparam.h param.h
BUILT_SOURCES = vzfft2d.f vzfft3d.f cparam.h

AM_CPPFLAGS = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a libfft_kernel.a $(LIBBENCH) @FLIBS@
doit.c: doit-ffte.c cparam.h

doitv_SOURCES=doitv.c
doitv_LDADD=libfftv.a libfft_kernel.a $(LIBBENCH) @FLIBS@
doitv.c: doit-ffte.c cparam.h

$(PKG)/param.h $(PKG)/fft235.f $(PKG)/mfft235.f $(PKG)/kernel.f $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f $(PKG)/vzfft1d.f $(PKG)/vzfft2d.f $(PKG)/vzfft3d.f: extract

cparam.h: $(PKG)/param.h extract
	rm -f $@
	sed 's,^C\(.*\)$$,/* \1 */,;s,[Pp][Aa][Rr][Aa][Mm][Ee][Tt][Ee][Rr] *(\([^)]*\)),static const int \1\;,' $(PKG)/param.h > $@

libfft_a_SOURCES = $(PKG)/factor.f $(PKG)/fft235.f $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f

libfftv_a_SOURCES = $(PKG)/factor.f $(PKG)/fft235.f $(PKG)/mfft235.f $(PKG)/vzfft1d.f vzfft2d.f vzfft3d.f

libfft_kernel_a_SOURCES = $(PKG)/kernel.f

zfft1d.o zfft2d.o zfft3d.o vzfft1d.o vzfft2d.o vzfft3d.o: param.h

# need a local copy or Fortran include statement is confused
param.h: $(PKG)/param.h
	rm -f $@
	$(CP) $(PKG)/param.h $@

# slightly patch to avoid F90 allocate statement, so we can compile with F77
vzfft2d.f: $(PKG)/vzfft2d.f vzfft2d.f.patch
	rm -f $@
	$(CP) $(PKG)/vzfft2d.f $@
	$(PATCH) $@ <$(srcdir)/vzfft2d.f.patch
vzfft3d.f: $(PKG)/vzfft3d.f vzfft3d.f.patch
	rm -f $@
	$(CP) $(PKG)/vzfft3d.f $@
	$(PATCH) $@ <$(srcdir)/vzfft3d.f.patch

include ../Makefile.common
