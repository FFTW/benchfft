if DOUBLE
PRG=doit
endif

PKG = numutils
TARFILE = numutil.tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH)

NUMUTILS = $(PKG)/commonlibs.c

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(NUMUTILS)
doit_LDADD=$(LIBBENCH)

$(NUMUTILS): extract

include ../Makefile.common
