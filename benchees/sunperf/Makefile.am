if HAVE_SUNPERF
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-lsunperf $(LIBBENCH) @FLIBS@

include ../Makefile.common

