if HAVE_MLIB
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) @MLIB_LIB@

include ../Makefile.common
