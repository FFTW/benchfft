if HAVE_MKL_DFTI
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-lmkl_ia32 $(LIBBENCH) @FLIBS@

include ../Makefile.common
