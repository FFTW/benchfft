if HAVE_SUNPERF
PRG=doit doit2d doit3d
endif

INCLUDES = $(INCLBENCH)

AM_LDFLAGS=-xlic_lib=sunperf

doit_SOURCES=doit.c common.h
doit_LDADD=$(LIBBENCH) @FLIBS@

doit2d_SOURCES=doit2d.c common.h
doit2d_LDADD=$(LIBBENCH) @FLIBS@

doit3d_SOURCES=doit3d.c common.h
doit3d_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common

