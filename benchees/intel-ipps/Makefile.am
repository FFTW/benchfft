if HAVE_IPPS
PRG=doit
endif

AM_CPPFLAGS = $(INCLBENCH)

doit_SOURCES=doit.c

doit_LDADD=-lipps -lippcore $(LIBBENCH)

include ../Makefile.common
