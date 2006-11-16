if HAVE_F77
if DOUBLE
if HAVE_SSE3
PRG=doit doitv doit_sse2 doitv_sse2 doit_sse3 doitv_sse3
noinst_LIBRARIES=libfft.a libfftv.a libfft_kernel.a libfft_sse2.a libfft_sse3.a
else
if HAVE_SSE2
PRG=doit doitv doit_sse2 doitv_sse2
noinst_LIBRARIES=libfft.a libfftv.a libfft_kernel.a libfft_sse2.a
else
PRG=doit doitv
noinst_LIBRARIES=libfft.a libfftv.a libfft_kernel.a
endif
endif
endif
endif

PKG = ffte-4.0
TARFILE = $(PKG).tgz
EXTRA_DIST = $(TARFILE) vzfft2d.f.patch vzfft3d.f.patch doit-ffte.c
LOCAL_COPIES = cparam.h param.h

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a libfft_kernel.a $(LIBBENCH) @FLIBS@
doit.c: doit-ffte.c cparam.h

doitv_SOURCES=doitv.c
doitv_LDADD=libfftv.a libfft_kernel.a $(LIBBENCH) @FLIBS@
doitv.c: doit-ffte.c cparam.h

doit_sse2_SOURCES=doit_sse2.c
doit_sse2_LDADD=libfft.a libfft_sse2.a $(LIBBENCH) @FLIBS@
doit_sse2.c: doit-ffte.c cparam.h

doitv_sse2_SOURCES=doitv_sse2.c
doitv_sse2_LDADD=libfftv.a libfft_sse2.a $(LIBBENCH) @FLIBS@
doitv_sse2.c: doit-ffte.c cparam.h

doit_sse3_SOURCES=doit_sse3.c
doit_sse3_LDADD=libfft.a libfft_sse3.a $(LIBBENCH) @FLIBS@
doit_sse3.c: doit-ffte.c cparam.h

doitv_sse3_SOURCES=doitv_sse3.c
doitv_sse3_LDADD=libfftv.a libfft_sse3.a $(LIBBENCH) @FLIBS@
doitv_sse3.c: doit-ffte.c cparam.h

$(PKG)/param.h $(PKG)/fft235.f $(PKG)/mfft235.f $(PKG)/sse2.c $(PKG)/sse3.c $(PKG)/kernel.f $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f $(PKG)/vzfft1d.f $(PKG)/vzfft2d.f $(PKG)/vzfft3d.f: extract

cparam.h: $(PKG)/param.h extract
	rm -f $@
	sed 's,^C\(.*\)$$,/* \1 */,;s,PARAMETER *(\([^)]*\)),static const int \1\;,' $(PKG)/param.h > $@

libfft_a_SOURCES = $(PKG)/fft235.f $(PKG)/zfft1d.f $(PKG)/zfft2d.f $(PKG)/zfft3d.f

libfftv_a_SOURCES = $(PKG)/fft235.f $(PKG)/mfft235.f $(PKG)/vzfft1d.f vzfft2d.f vzfft3d.f
BUILT_SOURCES = vzfft2d.f vzfft3d.f

libfft_kernel_a_SOURCES = $(PKG)/kernel.f

libfft_sse2_a_SOURCES = $(PKG)/sse2.c
libfft_sse2_a_CFLAGS = $(AM_CFLAGS) $(SSE2_CFLAGS)

libfft_sse3_a_SOURCES = $(PKG)/sse3.c
libfft_sse3_a_CFLAGS = $(AM_CFLAGS) $(SSE3_CFLAGS)

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
