if HAVE_ACML
PRG=doit
else
endif

INCLUDES = $(INCLBENCH)
BENCHMARK_USEROPT='-o=patient'

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) -lacml @FLIBS@

include ../Makefile.common
