if DOUBLE
PRG=doit
endif

PKG = mixfft05
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH)

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(PKG)/Mixfft.c
doit_LDADD=$(LIBBENCH) @FLIBS@

$(PKG)/Mixfft.c: extract

include ../Makefile.common
