if HAVE_MLIB
PRG=doit
endif

AM_CPPFLAGS = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) @MLIB_LIB@

include ../Makefile.common
