if SINGLE
if HAVE_F77
PRG=doit83 doit97 doit117
noinst_LIBRARIES=libfft.a
endif
endif

INCLUDES=$(INCLBENCH)

doit83_SOURCES=doit83.c
doit83_LDADD=libfft.a $(LIBBENCH) @FLIBS@

doit97_SOURCES=doit97.c
doit97_LDADD=libfft.a $(LIBBENCH) @FLIBS@

doit117_SOURCES=doit117.c
doit117_LDADD=libfft.a $(LIBBENCH) @FLIBS@

libfft_a_SOURCES=as117.f as83.f as97.f

include ../Makefile.common
