if SINGLE
if HAVE_F77
PRG=doit
endif
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c fft99f.f
doit_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common
