if HAVE_F77
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c fft.f
doit_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common
