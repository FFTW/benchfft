if HAVE_ESSL
PRG=doit
endif

AM_CPPFLAGS = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) -lessl $(FLIBS)

include ../Makefile.common
