if HAVE_SGIMATH
PRG=doit doit2d doit3d
else
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) @FLIBS@

doit2d_SOURCES=doit2d.c
doit2d_LDADD=$(LIBBENCH) @FLIBS@

doit3d_SOURCES=doit3d.c
doit3d_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common
