if DOUBLE
if HAVE_NAG
PRG = doit
endif
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES = doit.c
doit_LDADD = @LIBNAG@ $(LIBBENCH) @FLIBS@

include ../Makefile.common
