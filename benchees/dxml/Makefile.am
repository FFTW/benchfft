if HAVE_DXML
PRG=doit
else
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common
