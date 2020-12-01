PRG=doit
PKG=minfft-1.2.0
AM_CPPFLAGS = $(INCLBENCH) -I$(PKG)
TARFILE=$(PKG).tar.gz

if SINGLE
AM_CFLAGS=-DMINFFT_SINGLE
endif

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(PKG)/minfft.c
doit_LDADD=$(LIBBENCH)

doit.o: $(PKG)/minfft.h

$(PKG)/minfft.c $(PKG)/minfft.h: extract

include ../Makefile.common
