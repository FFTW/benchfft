if SINGLE
PRG=doit
endif

PKG=cwplib
INCLUDES = $(INCLBENCH) -I$(PKG)
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(PKG)/pfafft.c
doit_LDADD=$(LIBBENCH)

include ../Makefile.common

$(PKG)/pfafft.c $(PKG)/cwp.h: extract

doit.o: $(PKG)/cwp.h