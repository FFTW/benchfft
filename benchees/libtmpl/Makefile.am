if HAVE_LIBTMPL
if DOUBLE
PRG=doit
endif
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-ltmpl $(LIBBENCH)

include ../Makefile.common

