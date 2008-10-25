PRG=doit

INCLUDES = $(INCLBENCH) -I$(PKG)
PKG=FFTReal-2.00
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

doit_SOURCES=doit.cc
nodist_doit_SOURCES = $(PKG)/FFTReal.h
doit_LDADD=$(LIBBENCH)

include ../Makefile.common

$(PKG)/FFTReal.h: extract

doit.o: $(PKG)/FFTReal.h
