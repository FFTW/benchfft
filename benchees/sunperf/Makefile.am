if HAVE_SUNPERF
PRG=doit doit2d doit3d
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=-lsunperf $(LIBBENCH) @FLIBS@

doit2d_SOURCES=doit2d.c
doit2d_LDADD=-lsunperf $(LIBBENCH) @FLIBS@

doit3d_SOURCES=doit3d.c
doit3d_LDADD=-lsunperf $(LIBBENCH) @FLIBS@

include ../Makefile.common

