if DOUBLE
if HAVE_FFTW3
PRG=doit doitr
FFTW3_LIBS=@LFFTW3@
endif
else
if HAVE_FFTW3F
PRG=doit doitr
FFTW3_LIBS=@LFFTW3F@
endif
endif

INCLUDES = $(INCLBENCH)
BENCHMARK_USEROPT='-o=patient'

doit_SOURCES=doit.c
doit_LDADD=$(FFTW3_LIBS) $(LIBBENCH) @FLIBS@

doitr_SOURCES=doitr.c
doitr_LDADD=$(FFTW3_LIBS) $(LIBBENCH) @FLIBS@

include ../Makefile.common
