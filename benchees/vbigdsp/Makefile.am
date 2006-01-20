if SINGLE
PRG=doit
endif

PKG=vBigDSP
TARFILE=$(PKG)-20030322.tgz
EXTRA_DIST = $(TARFILE) vBigDSP.c.patch

INCLUDES = $(INCLBENCH) -I$(PKG)

AM_CFLAGS=@ALTIVEC_CFLAGS@ -Dvalloc=bench_malloc
doit_SOURCES=doit.c
nodist_doit_SOURCES = vBigDSP.c
doit_LDADD=$(LIBBENCH)

if INFO_ONLY
vBigDSP.c: extract
	rm -f $@
	echo "void FFTComplex(void){} void FFT2DComplex(void){} void FFTRealForward(void){} void FFT2DRealForward(void){} void FFTRealInverse(void){} void FFT2DRealInverse(void){} void ShutdownFFT(void) {}" > $@
else
vBigDSP.c: $(PKG)/vBigDSP.c
	rm -f $@
	$(CP) $(PKG)/vBigDSP.c $@
	$(PATCH) $@ <$(srcdir)/vBigDSP.c.patch

$(PKG)/vBigDSP.c: extract
endif
doit.o: extract

include ../Makefile.common
