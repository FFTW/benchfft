if SINGLE
if HAVE_F77
PRG=doit_sfftfu doit_ctfftsr
noinst_LIBRARIES=libfft.a
endif
endif

INCLUDES = $(INCLBENCH)

doit_sfftfu_SOURCES=doit_sfftfu.c 
doit_sfftfu_LDADD=libfft.a $(LIBBENCH) @FLIBS@

doit_ctfftsr_SOURCES=doit_ctfftsr.c 
doit_ctfftsr_LDADD=libfft.a $(LIBBENCH) @FLIBS@

libfft_a_SOURCES=csr1.f csr2.f icsrfft.f irvsrh.f rvsrh.f sfftfu.f

include ../Makefile.common
