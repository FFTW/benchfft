if SINGLE
PRG=doit
endif

PKG = ffts_in_C-1.0
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH)

KRUKAR = $(PKG)/b512.c $(PKG)/bitrev.c $(PKG)/dint.c $(PKG)/dintime.c	\
$(PKG)/idint.c $(PKG)/idintime.c $(PKG)/tab.c

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(KRUKAR)
doit_LDADD=$(LIBBENCH) @FLIBS@

$(KRUKAR): extract

include ../Makefile.common
