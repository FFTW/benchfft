if SINGLE
HARM=harm.f
EXTRA_DIST=harmd.f
else
HARM=harmd.f
EXTRA_DIST=harm.f
endif

if HAVE_F77
PRG=doit
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c $(HARM)
doit_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common
