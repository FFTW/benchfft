if HAVE_SGIMATH
PRG=doit doit2d doit3d
else
endif

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=$(LIBBENCH) -lcomplib.sgimath @FLIBS@

doit2d_SOURCES=doit2d.c
doit2d_LDADD=$(LIBBENCH) -lcomplib.sgimath @FLIBS@

doit3d_SOURCES=doit3d.c
doit3d_LDADD=$(LIBBENCH) -lcomplib.sgimath @FLIBS@

include ../Makefile.common
