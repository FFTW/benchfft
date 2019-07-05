if HAVE_VDSP
PRG=doit
endif

AM_CPPFLAGS = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) @VECLIB@

include ../Makefile.common
