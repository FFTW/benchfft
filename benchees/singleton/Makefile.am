if SINGLE
if HAVE_F77
PRG=doit doit_3d
noinst_LIBRARIES=libfft.a
endif
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@

doit_3d_SOURCES=doit_3d.c
doit_3d_LDADD=libfft.a $(LIBBENCH) @FLIBS@

libfft_a_SOURCES = fft.f realtr.f

include ../Makefile.common
