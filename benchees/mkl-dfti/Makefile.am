if HAVE_MKL_DFTI
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-lmkl $(LIBBENCH) @FLIBS@

include ../Makefile.common
