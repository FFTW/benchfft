if HAVE_VDSP
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) @VECLIB@

include ../Makefile.common
