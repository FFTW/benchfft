if HAVE_F77
if DOUBLE
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@

libfft_a_SOURCES = esrfft.f

include ../Makefile.common
