if HAVE_F77
PRG=doit doit_3d
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c fft.f
doit_LDADD=$(LIBBENCH) @FLIBS@

doit_3d_SOURCES=doit_3d.c fft.f
doit_3d_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common
