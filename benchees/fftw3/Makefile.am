if SINGLE
if HAVE_FFTW3F
PRG=doit
FFTW3_LIBS=@LFFTW3F@
endif
else
if HAVE_FFTW3
PRG=doit
FFTW3_LIBS=@LFFTW3@
endif
endif

INCLUDES = $(INCLBENCH)
BENCHMARK_USEROPT='-o=patient'

doit_SOURCES=doit.c
doit_LDADD=$(FFTW3_LIBS) $(LIBBENCH) @FLIBS@

include ../Makefile.common
