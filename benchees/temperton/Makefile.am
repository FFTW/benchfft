if SINGLE
if HAVE_F77
PRG=doit
noinst_LIBRARIES=libfft99f.a
endif
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c 
doit_LDADD=libfft99f.a $(LIBBENCH) @FLIBS@

libfft99f_a_SOURCES=fft99f.f

include ../Makefile.common
