if SINGLE
PRG=doit
endif

PKG=vBigDSP
TARFILE=$(PKG)-20030322.tgz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH) -I$(PKG)

CFLAGS=@CFLAGS@ @ALTIVEC_CFLAGS@ -Dvalloc=bench_malloc
doit_SOURCES=doit.c
nodist_doit_SOURCES = $(PKG)/vBigDSP.c
doit_LDADD=$(LIBBENCH)

if INFO_ONLY
$(PKG)/vBigDSP.c: extract
	rm -f $@
	echo "void FFTComplex(void){} void FFT2DComplex(void){} void FFTRealForward(void){} void FFT2DRealForward(void){} void FFTRealInverse(void){} void FFT2DRealInverse(void){} void ShutdownFFT(void) {}" > $@
else
$(PKG)/vBigDSP.c: extract
endif
doit.o: extract

include ../Makefile.common
