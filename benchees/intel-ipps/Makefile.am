if HAVE_IPPS
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-lipps $(LIBBENCH) @FLIBS@

include ../Makefile.common
