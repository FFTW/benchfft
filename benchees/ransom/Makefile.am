if SINGLE
PRG=doit
endif

PKG = ransom
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH) -I$(PKG)

RANSOM = $(PKG)/fastffts.c

doit_SOURCES=doit.c $(RANSOM)
BUILT_SOURCES=$(RANSOM)
doit_LDADD=$(LIBBENCH)

$(RANSOM): extract

include ../Makefile.common
