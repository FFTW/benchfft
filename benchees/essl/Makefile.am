if HAVE_ESSL
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) -lessl $(FLIBS)

include ../Makefile.common
